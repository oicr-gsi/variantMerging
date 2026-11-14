import argparse
import re
import traceback
import numpy as np

from pysam import bcftools as bt
from itertools import chain
from vcfParsing import *
from cyvcf2 import VCF, Writer

''' Check for sample swap (mutect2 - specific wrong order of NORMAL and TUMOR columns '''
def _samples_swapped(reader: cyvcf2.VCFReader, normalName: str, tumorName: str) -> bool:
    """ mutect2 - specific normal swap fix: """
    if len(reader.samples) == 2:
        # Mutect2 metadata-based detection
        if (
                (reader.samples[1] == 'NORMAL') or
                (
                        reader.contains('normal_sample') and
                        reader.samples[1] in reader['normal_sample']
                )
        ):
            return True

        # VarDict (or others): name-based detection
        elif normalName and tumorName:
            return reader.samples[0] == tumorName and reader.samples[1] == normalName
    return False

''' bcftools function for renaming only, re-header vcf and write into output file '''
def _rename_samples(input_path: str, output_path: str, normalName: str, tumorName: str):
    vcf_namefile = "sampleNames_fix.txt"
    with open(vcf_namefile, 'w') as f:
        f.write(f"NORMAL\t{normalName or 'NORMAL'}\n")
        f.write(f"TUMOR\t{tumorName or 'TUMOR'}\n")

    bt.reheader(
        '-s', vcf_namefile,
        '-o', output_path,
        input_path,
        catch_stdout=False
    )

''' Impute tumor/normal genotypes, GT value for records '''
def _tumor_normal_genotypes(ref, alt, info):
    """

    Retrieve standard 0/0, 0/1, 1/1 style genotypes from INFO field.
    Normal -- NT field (ref, het, hom, conflict)
    Tumor -- SGT field
      - for SNPs specified as GG->TT for the normal and tumor diploid alleles. These
        can also represent more complex alleles in which case we set at heterozygotes
        pending longer term inclusion of genotypes in Strelka2 directly
        (https://github.com/Illumina/strelka/issues/16)
      - For indels, uses the ref, het, hom convention
    ref: The REF allele from a VCF line
    alt: A list of potentially multiple ALT alleles (rec.ALT.split(";"))
    info: The VCF INFO field
    fname, coords: not currently used, for debugging purposes

    """
    known_names = {"het", "hom", "ref", "conflict"}

    def name_to_gt(val):
        if val.lower() == "het":
            return [0, 1, False]
        elif val.lower() == "hom":
            return [1, 1, False]
        elif val.lower() == "ref":
            return [0, 0, False]
        elif val.lower() == "conflict":
            return [-1, -1, False]
        else:
            """ Non-standard representations, het is our best imperfect representation """
            return [0, 1, False]

    def alleles_to_gt(val):
        gt_indices = {}
        for j, gT in enumerate([ref] + alt):
            gt_indices.update({str(gT).upper(): j})
        tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
        if tumor_gts and val not in known_names:
            if max(tumor_gts) == 0:
                t_gt = [0, 0, False]
            elif 0 in tumor_gts:
                t_gt = [0, min([x for x in tumor_gts if x > 0]), False]
            else:
                t_gt = [min(tumor_gts), max(tumor_gts), False]
        else:
            t_gt = name_to_gt(val)
        return t_gt

    nt_val = info['NT']
    normal_gt = name_to_gt(nt_val)
    sgt_val = info['SGT']
    if not sgt_val:
        tumor_gt = [0, 0, False]
    else:
        sgt_val = sgt_val.split("->")[-1]
        tumor_gt = alleles_to_gt(sgt_val)
    return [normal_gt, tumor_gt]


def _generate_ad(my_record: cyvcf2.Variant, samples: list):
    """
      Retrieve allele reads and generate AD in form of REF_count,ALT1_count,ALT2_count etc.
      ref: The REF allele from a VCF line
      alt: A list of potentially multiple ALT alleles (rec.ALT.split(";"))
    """
    ref = my_record.REF
    alt = my_record.ALT
    ref_count_ref = ref_count_alt = alt_count_ref = alt_count_alt = []
    try:
        if my_record.is_indel:
            ref_count_ref = [
                next(map(int, my_record.format('TAR')[0]))
            ]
            alt_count_ref = [
                next(map(int, my_record.format('TAR')[1]))
            ]
        else:
            ref_count_ref = [
                next(map(int, my_record.format(x + 'U')[0])) for x in ref
            ]

            alt_count_ref = [
                next(map(int, my_record.format(x + 'U')[1])) for x in ref
            ] if len(samples) > 1 else []

        if alt[0] is not None:
            if my_record.is_indel:
                ref_count_alt = [
                    next(map(int, my_record.format('TIR')[0]))
                ]
                alt_count_alt = [
                    next(map(int, my_record.format('TIR')[1]))
                ]
            else:
                ref_count_alt = [
                    next(map(int, my_record.format(x + 'U')[0])) for x in alt
                ]
                alt_count_alt = [
                    next(map(int, my_record.format(x + 'U')[1])) for x in alt
                ] if len(samples) > 1 else []
    except Exception:
        traceback.print_exc()
        print("Failed to extract AD value")
    ref_values = np.array([ref_count_ref, ref_count_alt], dtype=object)
    ref_values = list(chain.from_iterable(ref_values))
    alt_values = np.array([alt_count_ref, alt_count_alt], dtype=object)
    alt_values = list(chain.from_iterable(alt_values))
    return [ref_values, alt_values]


def main():
    parser = argparse.ArgumentParser(description='output file')
    parser.add_argument("path", type=str, help='input vcf file path')
    parser.add_argument('-o', '--output', help='output vcf file path', required=True)
    parser.add_argument('-r', '--reference', help='reference id, i.e. hg19', required=True)
    parser.add_argument('-t', '--tumor', help='tumor id, available for matched calls', required=False)
    parser.add_argument('-n', '--normal', help='normal id, available for matched calls', required=False)
    parser.add_argument('-f', '--renaming', action='store_true',
                        help='Fast mode - Only update sample names in header(no parsing of records)')

    args = parser.parse_args()

    """
     This script makes several things:

     * removes non-canonical contigs from both the header and records
     * adds GT and AD fields (imputed from NT, SGT and allele-specific read depths if available)
     * removes tool-specific header lines
     * changes sample names to NORMAL and TUMOR in header, reverse the order if TUMOR comes first
    """

    input_vcf = args.path
    output_vcf = args.output
    vcf_reader = VCF(input_vcf)

    """ 
       If NORMAL/TUMOR swap detected, correct it early
       Note that we cannot do swap correction and renaming in one go in some cases,
       so better do this in two passes if we have something unusual like
         * header shows sample names instead of NORMAL and TUMOR but the order is incorrect
       
       NORMAL/TUMOR in the column header allows both swap correction and renaming in one go
    """
    if _samples_swapped(vcf_reader, args.normal, args.tumor):
        vcf_swapfile = "./sampleSwap_fix.vcf.gz" if input_vcf.endswith('.gz') else "./sampleSwap_fix.vcf"
        samples_ordered = ",".join([vcf_reader.samples[1],vcf_reader.samples[0]])
        bt.view(
            '-s', samples_ordered,
            '-Oz' if input_vcf.endswith('gz') else '-Ov',
            '-o', vcf_swapfile,
            input_vcf,
            catch_stdout=False
        )
        vcf_reader.close()
        vcf_reader = VCF(vcf_swapfile)
        input_vcf = vcf_swapfile

    """
       Fast mode, we are only doing renaming if asked
       Note that this will cancel swap correction (won't be able to pass original normal and tumor names),
       we need to be sure that the order is correct when passing --renaming flag.
    """
    if args.renaming:
        print("We have sample names, will not parse record - only updating the header")
        vcf_reader.close()
        _rename_samples(input_vcf, output_vcf, args.normal, args.tumor)
        return

    """
     Process header lines:

     * Check for GT FORMAT, it needs to be present
     * get rid of alternative contigs, keep only canonical ones
     * use reference id passed as a parameter, do not use id derived from reference metadata field
     * Do not use any lines other than specified
    """
    add_gt = False
    have_filters = False

    ''' check if we have to add GT (required by downstream GATK analysis) '''
    try:
        add_gt = not any(
            (entry['HeaderType'] == 'FORMAT' and entry['ID'] == 'GT') for entry in vcf_reader.header_iter())
    except KeyError:
        print("Could not check if we have GT in FORMAT defined, no proper Key")
    ''' check if we have any filtering applied to records'''
    try:
        have_filters = any(entry['HeaderType'] == 'FILTER' for entry in vcf_reader.header_iter())
    except KeyError:
        print("Could not check if we have any FILTER defined, no proper Key")

    """
      Process Records:
      * Update header before opening writer
      * For records, make sure we have GT field properly formatted if we have SGT (strelka/strelka2-specific)
      * Only canonical contigs used in the pre-processed vcf
      * If record has no filters listed, we put 'PASS' if there were filters applied (or '.' if not)
    """

    update_header(vcf_reader, add_gt)

    """
       Open Writer, at this point we want all header updates complete
       Process all records
    """
    vcf_writer = Writer(args.output, vcf_reader)
    for record in vcf_reader:
        '''We process only canonical chromosomes'''
        if re.search('_', record.CHROM):
            continue

        '''Impute GT from available data (SGT)'''
        if add_gt and vcf_reader.contains('SGT'):
            gtStrings = _tumor_normal_genotypes(record.REF, '.' if record.ALT is None else record.ALT, record.INFO)
        if len(record.ALT) == 0 or record.ALT[0] is None:
            continue

        ''' ID, FILTER '''
        if record.ID is None:
            record.ID = '.'
        if record.FILTER is None or len(record.FILTER) == 0:
            record.FILTER = "." if not have_filters else 'PASS'

        ''' Adding GT and AD to a record '''
        if add_gt and gtStrings:
            gt_values = [
                gtStrings[0],
                gtStrings[1]
            ]
            ad_values = _generate_ad(record, vcf_reader.samples)
            record.genotypes = gt_values
            """ Ensure we have integer values for AD """
            ad_array = np.asarray(ad_values, dtype=np.int32)
            n_alleles = len(record.ALT) + 1
            assert ad_array.shape == (len(record.genotypes), n_alleles)
            record.set_format('AD', np.array(ad_array))
        vcf_writer.write_record(record)

    vcf_reader.close()
    vcf_writer.close()


if __name__ == '__main__':
    main()