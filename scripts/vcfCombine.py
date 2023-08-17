"""
Script for execution of two tasks -
* combine_vcfs
* consensus_vcf
"""
import argparse
import os
import VariantConsensus

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load list of vcf files to parse')
    parser.add_argument("list", type=str, help='list of vcf files to combine')
    parser.add_argument('-c', '--combined', help='combined vcf file path', required=False, default="combined.vcf")
    args = parser.parse_args()

    # Check if file exists
    if os.path.exists(args.list):
        print("Combining variants")
        vcf_list = {}
        for line in open(args.list, 'r'):
            vetted_line = line.strip()
            vcf_list[os.path.basename(vetted_line).rsplit(".", 2)[0]] = vetted_line
        variant_lines = VariantConsensus.combine_caller_vcfs(vcf_list)
        print("Writing into a file")
        with open(args.combined, 'w') as c:
            c.writelines(variant_lines)
