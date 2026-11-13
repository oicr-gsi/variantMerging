"""
   This module is a re-coded combine_vcfs.rb from VariantConsensus repo maintained by Miguel Vazquez (mikisvaz@gmail.com)
   check and correction of order of samples, pre-processing vcf header
   https://github.com/Rbbt-Workflows/VariantConsensus/blob/main/lib/VariantConsensus
"""
import gzip
import re
from cyvcf2 import VCF
from xopen import xopen

def _zip_fields(array: list, max_value=None):
    if array is None or len(array) == 0 or not any(array):
        return []

    first = array[0]
    if max_value is None:
        max_value = max(len(line) for line in array)

    rest = [v * max_value if v.length == 1 and max_value > 1 else v for v in array[1:]]
    if len(first) == 1 and max_value > 1:
        first = first * max_value

    return list(zip(first, *rest))


def zip_fields(array: list):
    if len(array) < 10000:
        return _zip_fields(array)
    else:
        zipped_slices = []
        max_value = max(len(line) for line in array)
        for slice_ in (array[i:i+10000] for i in range(0, len(array), 10000)):
            zipped_slices.append(_zip_fields(slice_, max_value))

        new = zipped_slices[0]
        for rest in zipped_slices[1:]:
            for i, list_ in enumerate(rest):
                new[i].extend(list_)

        return new


def get_header_lines(file_path: str):
    vcf_reader = VCF(file_path)
    header_lines = vcf_reader.raw_header
    vcf_reader.close()
    lines = header_lines.split("\n")
    return list(filter(None, map(str.strip, lines)))


def guess_vcf_tumor_sample(vcf: str):
    header_lines = get_header_lines(vcf)
    for line in header_lines:
        match = re.search(r"tumor_sample=(\S+)", line)
        if match:
            return match.group(1)
    fields = header_lines[-1].split("\t")
    sample1 = fields[-2]
    if len(fields) == 9:
        print(f"Could not find tumor_sample field in {vcf}, but only one sample: {fields[-1]}")
        return fields[-1]
    elif "TUMOR" in fields or "Tumor" in fields:
        print(f"Could not find tumor_sample field in {vcf}, using TUMOR")
        return "TUMOR"
    else:
        try:
            vcf_reader = VCF(vcf)
            for entry in vcf_reader:
                if entry and any(
                        genotype in ["0/1", "0|1", "1/1", "1|1"] for genotype in entry[sample1].split(":")):
                    genotype = next(genotype for genotype in entry[sample1].split(":") if
                                    genotype in ["0/1", "0|1", "1/1", "1|1"])
                    print(f"Could not find tumor_sample field in {vcf}, but {sample1} has genotype {genotype}")
                    vcf_reader.close()
                    return sample1
            vcf_reader.close()
        except:
            print("error parsing genotype")
        print(f"Could not find tumor_sample field in {vcf}, using last field")
        return fields[-1]


def combine_caller_vcfs(call_list: dict, output_file: str):
    preambles = {}
    for name, file in call_list.items():
        preambles[name] = get_header_lines(file)

    preamble = []
    fields = None
    headers = {}

    for name, lines in preambles.items():
        for line in lines:
            if line.startswith("#CHR"):
                headers[name] = line.split("\t")
                continue

            line = line.replace("FILTER=<ID=", f"FILTER=<ID={name}--")
            line = line.replace("FORMAT=<ID=", f"FORMAT=<ID={name}--")
            line = line.replace("INFO=<ID=", f"INFO=<ID={name}--")
            if line not in preamble and not any(keyword in line for keyword in ["tumor_sample", "normal_sample"]):
                preamble.append(line)

    preamble.append('##INFO=<ID=CalledBy,Number=.,Type=String,Description="Callers calling this variant">')

    for name in call_list.keys():
        if not any(line for line in preamble if "FILTER" in line and f"{name}--PASS" in line):
            preamble.append(f"##FILTER=<ID={name}--PASS,Description=\"Passes {name} filter\">")

    preamble.extend([
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position in the sample">',
        '##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">',
        '##FORMAT=<ID=AF,Number=1,Type=String,Description="Allele fractions of alternate alleles in the tumor">',
        '##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of A alleles used in tiers 1,2">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Depth of reads supporting alleles or variant allele">',
        '##FORMAT=<ID=BCOUNT,Number=4,Type=Integer,Description="Occurrence count for each base at this site (A,C,G,T)">'
    ])
    variants = {}
    called_by = {}
    for name, file in call_list.items():
        tumor_id = guess_vcf_tumor_sample(file)
        swap_samples = tumor_id == headers[name][-2] or headers[name][-1] == "NORMAL"
        fields = headers[name]
        fields[-1] = "TUMOR"
        fields[-2] = "NORMAL"
        with xopen(file, "rt", encoding="utf-8") as vr:
            for line in vr:
                if line.startswith("#"):
                    continue
                line = line.rstrip()
                # TODO: introduce a support for tumor-only data in our script sometime in a future
                chrom, pos, rsid, ref, alt, qual, vfilter, info, vformat, normal_sample, tumor_sample, *rest = line.split('\t')
                if swap_samples:
                    normal_sample, tumor_sample = tumor_sample, normal_sample
                mutation = ":".join([chrom, pos, ref, alt])
                vfilter = ";".join([f"{name}--{f}" for f in vfilter.split(";") if f != "."])
                info = ";".join([f"{name}--{f}" for f in info.split(";") if f != "."])
                vformat = ":".join([f"{name}--{f}" for f in vformat.split(":") if f != "."])
                variants.setdefault(mutation, []).append([vfilter, info, vformat, normal_sample, tumor_sample] + rest)
                called_by.setdefault(mutation, []).append(name)

    str_ = "\n".join(preamble) + "\n"
    str_ += "\t".join(fields) + "\n"

    for mutation, lists in variants.items():
        chrom, pos, ref, alt = mutation.split(":")
        mfilter = []
        minfo = []
        mformat = []
        msample1 = []
        msample2 = []
        mrest = []
        for entry in lists:
            filter_, info, vformat, sample1, sample2, *rest = entry
            mfilter.append(filter_)
            minfo.append(info)
            mformat.extend(vformat.split(":"))

            if sample2 is None:
                msample2.extend(sample1.split(":"))
                msample1.extend(["."] * len(sample1.split(":")))
            else:
                msample1.extend(sample1.split(":"))
                msample2.extend(sample2.split(":"))
            mrest.append(rest)

        new_format = []
        new_sample1 = []
        new_sample2 = []

        common_format = ["GT", "AF", "DP", "AU", "AD", "BCOUNT", "DP4", "TIR"]
        for key in common_format:
            match = next((k for k in mformat if k.split("--")[-1] == key), None)
            if match:
                kpos = mformat.index(match)
                new_format.append(match.partition("--")[-1])
                new_sample1.append(msample1[kpos])
                new_sample2.append(msample2[kpos])

        mformat = new_format + mformat
        msample1 = new_sample1 + msample1
        msample2 = new_sample2 + msample2

        minfo = [i for i in minfo if i]
        minfo.insert(0, f"CalledBy={','.join(called_by[mutation])}")
        str_ += "\t".join([chrom, pos, ".", ref, alt, ".", ";".join(mfilter), ";".join(minfo), ":".join(mformat),
                           ":".join(msample1), ":".join(msample2)] + [":".join(r) for r in zip_fields(mrest)]) + "\n"

    if output_file.endswith(".gz"):
        with gzip.open(output_file, "wb") as f:
            f.write(str_.encode("utf-8"))
    else:
        with open(output_file, "+w") as f:
            f.writelines(str_)
