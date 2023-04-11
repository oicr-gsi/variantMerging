"""
   Module for operations on vcf files using PyVCF package
   the goal is to reduce redundancy  - move code used by scripts for
   pre- and post-processing into one place

"""

import vcf
import re

def _special_character(f):
    switcher = {-1: 'A',
                -2: 'G',
                -3: 'R'}
    return switcher[f] if f in switcher.keys() else '.'


def generate_header(vcf_reader: vcf.Reader, reference: str):
    header_lines = ["##fileformat=" + vcf_reader.metadata['fileformat'] + "\n"]
    add_gt = True

    for flt in vcf_reader.filters:
        header_lines.append("##FILTER=<ID=" + vcf_reader.filters[flt].id +
                            ",Description=\"" + vcf_reader.filters[flt].desc + "\">\n")
    '''
       Process FORMATs
    '''
    for f in vcf_reader.formats:
        num = vcf_reader.formats[f].num
        if vcf_reader.formats[f].num is None:
            num = '.'
        elif isinstance(vcf_reader.formats[f].num, int) and vcf_reader.formats[f].num < 0:
            num = _special_character(vcf_reader.formats[f].num)
        if vcf_reader.formats[f].id == 'GT':
            add_gt = False
        header_lines.append("##FORMAT=<ID=" + vcf_reader.formats[f].id +
                            ",Number=" + str(num) + ",Type=" + vcf_reader.formats[f].type +
                            ",Description=\"" + vcf_reader.formats[f].desc + "\">\n")
    if add_gt:
        header_lines.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\""
                            "Genotype, constructed from SGT INFO via external modification\">\n")
        header_lines.append("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\""
                            "Allelic depths for the ref and alt alleles in the order listed\">\n")

    ''' 
       Process INFO
    '''
    for i in vcf_reader.infos:
        num = vcf_reader.infos[i].num
        if vcf_reader.infos[i].num is None:
            num = '.'
        elif isinstance(vcf_reader.infos[i].num, int) and vcf_reader.infos[i].num < 0:
            num = _special_character(vcf_reader.infos[i].num)
        header_lines.append("##INFO=<ID=" + vcf_reader.infos[i].id +
                            ",Number=" + str(num) + ",Type=" + vcf_reader.infos[i].type +
                            ",Description=\"" + vcf_reader.infos[i].desc + "\">\n")

    ''' 
        Process contigs
    '''

    for c in vcf_reader.contigs:
        if re.search('_', vcf_reader.contigs[c].id):
            continue
        header_lines.append("##contig=<ID=" +
                            vcf_reader.contigs[c].id +
                            ",length=" +
                            str(vcf_reader.contigs[c].length) +
                            ",assembly=" + reference + ">\n")

    header_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']
    header_lines.append("#" + "\t".join(header_fields) + "\n")
    return header_lines



