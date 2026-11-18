"""
   Module for operations on vcf files using cyvcf2 package
   https://github.com/brentp/cyvcf2
   the goal is to reduce redundancy  - move code used by scripts for
   pre- and post-processing into one place

"""
import cyvcf2

def _special_character(f):
    switcher = {-1: 'A',
                -2: 'G',
                -3: 'R'}
    return switcher[f] if f in switcher.keys() else '.'

def update_header(vcf_reader: cyvcf2.Reader, add_gt: bool):
    if add_gt:
        vcf_reader.add_format_to_header({
           'ID': 'GT',
           'Description': 'Genotype, constructed from SGT INFO via external modification',
           'Type': 'Character',
           'Number': '1'
        })
        vcf_reader.add_format_to_header({
            'ID': 'AD',
            'Description': 'Allelic depths for the ref and alt alleles in the order listed',
            'Type': 'Integer',
            'Number': 'R'
        })
    ''' Handle negative Number instances '''
    for field in vcf_reader.header_iter():
        """ Process FORMAT and INFO"""
        if field['HeaderType'] in ['FORMAT', 'INFO']:
            if not field['Number']:
                vcf_reader.update(field['ID'], field['Type'], '.', field['Description'])
            elif isinstance(field['Number'], int) and field['Number'] < 0:
                vcf_reader.update(field['ID'], field['Type'], _special_character(field['NUM']), field['Description'])

