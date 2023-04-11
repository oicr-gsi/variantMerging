import argparse
from vcfParsing import *

parser = argparse.ArgumentParser(description='output file')
parser.add_argument("-m", '--merged', type=str, help='merged vcf file path', required=True)
parser.add_argument("-c", '--combined', type=str, help='combined vcf file path', required=True)
parser.add_argument("-r", "--reference", type=str, help='reference id, such as hg38', required=True)
parser.add_argument('-o', '--output', type=str, help='output vcf file path', required=True)
args = parser.parse_args()

'''
   A script which corrects the 'set' tag according to the overlap of variants
   from different callers. DISCVRseq currently uses names of workflows but more complex overlaps
   are not annotated properly
'''

merged_reader = vcf.Reader(filename=args.merged, compressed=True)
combined_reader = vcf.Reader(filename=args.combined, compressed=True)

'''
  Use filter field info to construct tags for call made by more than one caller:

  1. intersect == PASS in all callers
  2. filtered_in_all == no PASS calls in any of the callers
  3. partially_filtered == some are PASS, some filtered
'''

def _determine_set(quality_hash: list):
    if "" in quality_hash['filters']:
        for q in quality_hash['filters']:
            if q is not "":
                return "partially_filtered"
        return "intersect"
    else:
        return "filtered_in_all"

record_lines = []
already_seen = {}

for record in merged_reader:
    my_id = "_".join([record.CHROM, str(record.POS), "".join(map(str, record.alleles))])
    if my_id in already_seen.keys():
        already_seen[my_id]['filters'].append(":".join(map(str, record.FILTER)))
        already_seen[my_id]['count'] += 1
    else:
        already_seen[my_id] = {'count': 1, 'filters': [":".join(map(str, record.FILTER))]}

'''
 Process header lines:
 * get rid of alternative contigs, keep only canonical ones
 * Do not use any lines other than specified
'''

header_lines = generate_header(combined_reader, args.reference)

'''
   Process records in combined (by DISCVRseq) vcf, change set appropriately if needed
'''
for record in combined_reader:
    my_id = "_".join([record.CHROM, str(record.POS), "".join(map(str, record.alleles))])
    new_set = None
    if my_id in already_seen.keys() and already_seen[my_id]['count'] > 1:
        new_set = _determine_set(already_seen[my_id])

    record_data = [record.CHROM, str(record.POS)]
    idString = "." if record.ID is None else record.ID
    qString = "." if record.QUAL is None else record.QUAL
    altString = ",".join(map(str, record.ALT))
    if record.FILTER is None or len(record.FILTER) == 0:
        filterString = "." if len(combined_reader.filters) == 0 else 'PASS'
    else:
        filterString = ";".join(record.FILTER)

    record_data.extend([idString, record.REF, altString, qString, filterString])

    info_data = []
    for field in record.INFO:
        if isinstance(record.INFO[field], list):
            info_data.append("=".join([field, ",".join(map(str, record.INFO[field]))]))
        elif field in combined_reader.infos.keys():
            f_type = combined_reader.infos[field].type
            if f_type == 'Flag' and record.INFO[field]:
                info_data.append(field)
            elif field == "set" and new_set:
                info_data.append("=".join([field, new_set]))
            else:
                info_data.append("=".join([field, "." if record.INFO[field] is None else str(record.INFO[field])]))

    info_string = ";".join(info_data) if len(info_data) > 0 else "."
    record_data.append(info_string)

    record_lines.append("\t".join(record_data) + "\n")


with open(args.output, mode='+w') as out:
    out.writelines(header_lines)
    out.writelines(record_lines)
out.close()
