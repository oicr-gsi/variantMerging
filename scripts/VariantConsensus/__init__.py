"""
This module is a re-coded combine_vcfs.rb from VariantConsensus repo maintained by Miguel Vazquez (mikisvaz@gmail.com)
check and correction of order of samples, pre-processing vcf header
"""
import csv
import gzip
import subprocess


def _zip_fields(array, max_value=None):
    # TODO: ask Miguel about this (checking nested array for emptiness?)
    if array is None or len(array) == 0 or not any(array):
        return []

    first = array[0]
    if max_value is None:
        max_value = max(len(l) for l in array)

    rest = [v * max_value if v.length == 1 and max_value > 1 else v for v in array[1:]]
    if len(first) == 1 and max_value > 1:
        first = first * max_value

    return list(zip(first, *rest))


def zip_fields(array):
    if len(array) < 10000:
        return _zip_fields(array)
    else:
        zipped_slices = []
        max_value = max(len(l) for l in array)
        for slice_ in (array[i:i+10000] for i in range(0, len(array), 10000)):
            zipped_slices.append(_zip_fields(slice_, max_value))

        new = zipped_slices[0]
        for rest in zipped_slices[1:]:
            for i, list_ in enumerate(rest):
                new[i].extend(list_)

        return new


def get_comment_lines(file_path):
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as file:
            comment_lines = [line.strip() for line in file if line.startswith('#')]
            return comment_lines
    else:
        with open(file_path, 'r') as file:
            comment_lines = [line.strip() for line in file if line.startswith('#')]
            return comment_lines


def get_data_lines(file_path):
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as file:
            data_lines = [line.strip() for line in file if not line.startswith('#')]
            return data_lines
    else:
        with open(file_path, 'r') as file:
            data_lines = [line.strip() for line in file if not line.startswith('#')]
            return data_lines

def guess_vcf_tumor_sample(vcf):
    try:
        if vcf.endswith('.gz'):
            grep_output = subprocess.check_output(f"zgrep 'tumor_sample=' '{vcf}'", shell=True).decode().strip()
        else:
            grep_output = subprocess.check_output(f"grep 'tumor_sample=' '{vcf}'", shell=True).decode().strip()
        return grep_output.split("=")[-1]
    except subprocess.CalledProcessError:
        header_lines = get_comment_lines(vcf)
        fields = header_lines[-1].split("\t")
        sample1, sample2 = fields[-2], fields[-1]
        if len(fields) == 9:
            print(f"Could not find tumor_sample field in {vcf}, but only one sample: {fields[-1]}")
            return fields[-1]
        elif "TUMOR" in fields or "Tumor" in fields:
            print(f"Could not find tumor_sample field in {vcf}, using TUMOR")
            return "TUMOR"
        else:
            data_lines = get_data_lines(vcf)
            for entry in data_lines:
                if entry and any(
                        genotype in ["0/1", "0|1", "1/1", "1|1"] for genotype in entry[sample1].split(":")):
                    genotype = next(genotype for genotype in entry[sample1].split(":") if
                                    genotype in ["0/1", "0|1", "1/1", "1|1"])
                    print(f"Could not find tumor_sample field in {vcf}, but {sample1} has genotype {genotype}")
                    return sample1
            print(f"Could not find tumor_sample field in {vcf}, using last field")
            return fields[-1]


def combine_caller_vcfs(list):
    preambles = {}
    for name, file in list.items():
        preambles[name] = get_comment_lines(file)
        #with open(file) as f:
        #    preambles[name] = [line for line in f.read().split("\n") if line.startswith("#")]

    preamble = []

    fields = None

    for name, lines in preambles.items():
        for line in lines:
            if line.startswith("#CHR"):
                fields = line.split("\t")
                continue

            line = line.replace("FILTER=<ID=", f"FILTER=<ID={name}--")
            line = line.replace("FORMAT=<ID=", f"FORMAT=<ID={name}--")
            line = line.replace("INFO=<ID=", f"INFO=<ID={name}--")
            if line not in preamble and not any(keyword in line for keyword in ["tumor_sample", "normal_sample"]):
                preamble.append(line)

    preamble.append('##INFO=<ID=CalledBy,Number=.,Type=String,Description="Callers calling this variant">')

    for name in list.keys():
        if not any(l for l in preamble if "FILTER" in l and f"{name}--PASS" in l):
            preamble.append(f"##FILTER=<ID={name}--PASS,Description=\"Passes {name} filter\">")

    preamble.extend([
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position in the sample">',
        '##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">',
        '##FORMAT=<ID=AF,Number=1,Type=String,Description="Allele fractions of alternate alleles in the tumor">',
        '##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of A alleles used in tiers 1,2">',
        '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depth of reads supporting alleles or variant allele">',
        '##FORMAT=<ID=BCOUNT,Number=4,Type=Integer,Description="Occurrence count for each base at this site (A,C,G,T)">'
    ])

    variants = {}
    called_by = {}
    for name, file in list.items():
        # TODO ask Miguel about this (why we are getting this if tumor_sample gets set below)
        tumor_id = guess_vcf_tumor_sample(file)
        vcf_lines = get_data_lines(file)
        for line in vcf_lines:
            # TODO: ask Miguel about this (is the intention to process matched calls only?)
            chrom, pos, rsid, ref, alt, qual, vfilter, info, vformat, normal_sample, tumor_sample, *rest = line.split('\t')
            swap_samples = tumor_id == fields[-2] or fields[-1] == "NORMAL"
            # TODO : ask Miguel about this (how do we check for swapped TUMOR and NORMAL)
            #if len(vformat.split(":")) == 10:
            #    sample1, sample2 = vformat.split(":")[-2:]
            #    swap_samples: bool = sample1 == tumor_sample
            #else:
            #    sample2 = vformat.split(":")[-1]
            #    swap_samples = False
            if swap_samples:
                normal_sample, tumor_sample = tumor_sample, normal_sample
            mutation = ":".join([chrom, pos, ref, alt])
            vfilter = ";".join([f"{name}--{f}" for f in vfilter.split(";") if f != "."])
            info = ";".join([f"{name}--{f}" for f in info.split(";") if f != "."])
            vformat = ":".join([f"{name}--{f}" for f in vformat.split(":") if f != "."])
            variants.setdefault(mutation, []).append([vfilter, info, vformat, normal_sample, tumor_sample] + rest)
            called_by.setdefault(mutation, []).append(name)

    fields[-2] = "NORMAL"
    fields[-1] = "TUMOR"

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

    return str_
