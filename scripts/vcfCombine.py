"""
   Script organizes inputs and execute VariantConsensus task -
   * combine_caller_vcfs
"""
import argparse
import os
import VariantConsensus

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load list of vcf files to parse')
    parser.add_argument("list", type=str, help='list of vcf files to combine')
    parser.add_argument("-n", "--names", help="Names of the callers, same order as in the list", required=True)
    parser.add_argument('-c', '--combined', help='combined vcf file path', required=False, default="combined.vcf")
    args = parser.parse_args()

    if os.path.exists(args.list):
        print("Combining variants")
        vcf_list = {}
        names_index = 0
        wf_names = args.names.split(",")
        for line in open(args.list, 'r'):
            vetted_line = line.strip()
            vcf_list[wf_names[names_index]] = vetted_line
            names_index += 1
        variant_lines = VariantConsensus.combine_caller_vcfs(vcf_list)
        print("Writing into a file")
        with open(args.combined, 'w') as c:
            c.writelines(variant_lines)
