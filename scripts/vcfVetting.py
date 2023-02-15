import vcf
import re
import argparse

parser = argparse.ArgumentParser(description='output file')
parser.add_argument("path", type=str, help='input vcf file path')
parser.add_argument('-o', '--output', help='output vcf file path', required=True)
parser.add_argument('-r', '--reference', help='reference id, i.e. hg19', required=True)
args = parser.parse_args()

'''
 This script makes several things:

 * removes non-canonical contigs from both the header and records
 * adds GT and AD fields (imputed from NT, SGT and allele-specific read depths if available)
 * removes tool-specific header lines
 * changes sample names to NORMAL and TUMOR in header, reverse the order if TUMOR comes first
'''

vcf_reader = vcf.Reader(filename=args.path, compressed=True)
header_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
header_lines = ["##fileformat=" + vcf_reader.metadata['fileformat'] + "\n"]

'''
 vcf 4.2 specs: special character handling
'''

def _special_character(f):
    switcher = {-1: 'A',
                -2: 'G',
                -3: 'R'}
    return switcher[f] if f in switcher.keys() else '.'


'''
 Process header lines:

 * Check for GT FORMAT, it needs to be present
 * get rid of alternative contigs, keep only canonical ones
 * use reference id passed as a parameter, do not use id derived from reference metadata field
 * Do not use any lines other than specified
'''

for f in vcf_reader.filters:
    header_lines.append("##FILTER=<ID=" + vcf_reader.filters[f].id +
                        ",Description=\"" + vcf_reader.filters[f].desc + "\">\n")
add_gt = True    # Add GT and AD formats if missing (strelka-specific)
swap_nt = False  # Swap sample data if we have TUMOR preceding NORMAL (mutect2-specific)

for m in vcf_reader.formats:
    num = vcf_reader.formats[m].num
    if vcf_reader.formats[m].num is None:
        num = '.'
    elif vcf_reader.formats[m].num < 0:
        num = _special_character(vcf_reader.formats[m].num)
    if vcf_reader.formats[m].id == 'GT':
        add_gt = False
    header_lines.append("##FORMAT=<ID=" + vcf_reader.formats[m].id +
                        ",Number=" + str(num) + ",Type=" + vcf_reader.formats[m].type +
                        ",Description=\"" + vcf_reader.formats[m].desc + "\">\n")
if add_gt:
    header_lines.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\""
                        "Genotype, constructed from SGT INFO via external modification\">\n")
    header_lines.append("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\""
                        "Allelic depths for the ref and alt alleles in the order listed\">\n")

for i in vcf_reader.infos:
    num = vcf_reader.infos[i].num
    if vcf_reader.infos[i].num is None:
        num = '.'
    header_lines.append("##INFO=<ID=" + vcf_reader.infos[i].id +
                        ",Number=" + str(num) + ",Type=" + vcf_reader.infos[i].type +
                        ",Description=\"" + vcf_reader.infos[i].desc + "\">\n")

for c in vcf_reader.contigs:
    if re.search('_', vcf_reader.contigs[c].id):
        continue
    header_lines.append("##contig=<ID=" + vcf_reader.contigs[c].id + ",length=" + str(vcf_reader.contigs[c].length) +
                        ",assembly=" + args.reference + ">\n")

inputHash = {}

'''
    Preserve metadata for inputs
'''
if 'inputs' in vcf_reader.metadata.keys():
    header_lines.append("##inputs=" + " ".join(vcf_reader.metadata['inputs']) + "\n")
    for s in vcf_reader.metadata['inputs'][0].split(" "):
        keyValuePair = s.split(":")
        inputHash.update({keyValuePair[1]: keyValuePair[0]})
    sampleList = list(inputHash.values())
    if sampleList[1] == 'NORMAL':
        swap_nt = True

""" mutect2 - specific normal swap fix: """
if len(vcf_reader.samples) == 2 and vcf_reader.samples[1] == 'NORMAL' and len(inputHash) == 0:
    swap_nt = True
if 'normal_sample' in vcf_reader.metadata.keys() and vcf_reader.samples[1] in vcf_reader.metadata['normal_sample']:
    swap_nt = True

if vcf_reader.samples[0] not in ['NORMAL', 'TUMOR'] and len(inputHash) == 2:
    header_fields.extend(inputHash.values() if not swap_nt else ['NORMAL', 'TUMOR'])
else:
    header_fields.extend(vcf_reader.samples if not swap_nt else reversed(vcf_reader.samples))
header_lines.append("#" + "\t".join(header_fields) + "\n")

'''
  Definitions for helper subroutines which 

  * Process NT/SGT fields from strelka an generate GT
  * Generate AD numbers (also with strelka)
'''

def _sgtToGt(formatString, ref, alt):
    sgtStrings = formatString['SGT'].split("->")
    gt_norm = '0/0' if sgtStrings[0] == ref + ref else '1/1' if sgtStrings[0] == alt + alt else '0/1'
    gt_tumor = '0/0' if sgtStrings[1] == ref + ref else '1/1' if sgtStrings[1] == alt + alt else '0/1'
    return [gt_norm, gt_tumor]

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
    known_names = set(["het", "hom", "ref", "conflict"])

    def name_to_gt(val):
        if val.lower() == "het":
            return "0/1"
        elif val.lower() == "hom":
            return "1/1"
        elif val.lower() in set(["ref", "conflict"]):
            return "0/0"
        else:
            """ Non-standard representations, het is our best imperfect representation """
            return "0/1"

    def alleles_to_gt(val):
        gt_indices = {}
        for j, gT in enumerate([ref] + alt):
            gt_indices.update({str(gT).upper(): j})
        tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
        if tumor_gts and val not in known_names:
            if max(tumor_gts) == 0:
                tumor_gt = "0/0"
            elif 0 in tumor_gts:
                tumor_gt = "0/%s" % min([x for x in tumor_gts if x > 0])
            else:
                tumor_gt = "%s/%s" % (min(tumor_gts), max(tumor_gts))
        else:
            tumor_gt = name_to_gt(val)
        return tumor_gt

    nt_val = info['NT']
    normal_gt = name_to_gt(nt_val)
    sgt_val = info['SGT']
    if not sgt_val:
        tumor_gt = "0/0"
    else:
        sgt_val = sgt_val.split("->")[-1]
        tumor_gt = alleles_to_gt(sgt_val)
    return tumor_gt, normal_gt


def _generate_ad(ref, alt, samples):
    """
      Retrieve allele reads and generate AD in form of REF_count,ALT1_count,ALT2_count etc.
      ref: The REF allele from a VCF line
      alt: A list of potentially multiple ALT alleles (rec.ALT.split(";"))
    """
    ref_count_ref = ','.join([str(samples[0][x + 'U'][0]) for x in ref])
    alt_count_ref = ','.join([str(samples[1][x + 'U'][0]) for x in ref])
    ref_count_alt = "0"
    alt_count_alt = "0"
    if alt[0] is not None:
        ref_count_alt = ','.join([str(samples[0][str(x) + 'U'][0]) for x in alt])
        alt_count_alt = ','.join([str(samples[1][str(x) + 'U'][0]) for x in alt])
    ref_values = [ref_count_ref, ref_count_alt]
    alt_values = [alt_count_ref, alt_count_alt]
    return_values = [ref_values, alt_values] if record.samples[0].sample == 'NORMAL' else [alt_values, ref_values]
    return return_values


'''
  Process Records:
  * For records, make sure we have GT field properly formatted if we have SGT (strelka/strelka2-specific)
  * Only canonical contigs used in the pre-processed vcf
  * If record has no filters listed, we put 'PASS' if there were filters applied (or '.' if not)
'''
record_lines = []
for record in vcf_reader:
    if re.search('_', record.CHROM):
        continue
    if add_gt and 'SGT' in vcf_reader.infos.keys():
        gtStrings = _tumor_normal_genotypes(record.REF, '.' if record.ALT is None else record.ALT, record.INFO)
    if record.ALT[0] is None:
        continue

    record_data = [record.CHROM, str(record.POS)]
    idString = "." if record.ID is None else record.ID
    qString = "." if record.QUAL is None else record.QUAL
    altString = ",".join(map(str, record.ALT))
    if record.FILTER is None or len(record.FILTER) == 0:
        filtString = "." if len(vcf_reader.filters) == 0 else 'PASS'
    else:
        filtString = ";".join(record.FILTER)
    record_data.extend([idString, record.REF, altString, qString, filtString])

    """ Process INFO values, Flag type needs to be printed without value """
    info_data = []
    for field in record.INFO:
        if isinstance(record.INFO[field], list):
            info_data.append("=".join([field, ",".join(map(str, record.INFO[field]))]))
        elif field in vcf_reader.infos.keys():
            f_type = vcf_reader.infos[field].type
            if f_type == 'Flag' and record.INFO[field]:
                info_data.append(field)
            else:
                info_data.append("=".join([field, "." if record.INFO[field] is None else str(record.INFO[field])]))
    info_string = ";".join(info_data) if len(info_data) > 0 else "."
    record_data.append(info_string)

    """ Add FORMAT keys """
    format_data = ['GT:AD'] if add_gt else []
    for field in record.FORMAT.split(":"):
        format_data.append(field)
    record_data.append(":".join(format_data))

    """ Process FORMAT values """
    if add_gt:
        gt_values = [gtStrings[1], gtStrings[0]] if record.samples[0].sample == 'NORMAL' else gtStrings
        ad_values = _generate_ad([record.REF], record.ALT, record.samples)
    sample_data = []
    for sample in record.samples:
        format_data = [gt_values.pop(0), ",".join(ad_values.pop(0))] if add_gt else []
        for field in record.FORMAT.split(":"):
            if isinstance(sample[field], list):
                format_data.append(",".join(map(str, sample[field])))
            else:
                format_data.append("." if sample[field] is None else str(sample[field]))
        sample_data.append(":".join(format_data))
    if swap_nt:
        sample_data.reverse()
    for s in sample_data:
        record_data.append(s)
    record_lines.append("\t".join(record_data) + "\n")

'''
   PRINTING OUTPUT: the output file name is passed via -o option
'''

with open(args.output, mode='+w') as out:
    out.writelines(header_lines)
    out.writelines(record_lines)
out.close()
