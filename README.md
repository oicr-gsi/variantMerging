# variantMerging

VariantMerging is a workflow for combining variant calls from SNV analyses done with different callers (such as muTect2, strelka2). The workflow pre-processes input vcf files by removing non-canonical contigs, fixing fields and inferring missing values from available data. It combines calls, annotating them with caller-specific tags which allows identification of consensus variants. The workflow also uses GATK for producing merged results. In this case, all calls appear as-as. Essentially, this is a simple concatenation of the inputs.

### Pre-processing

The script used at this step performs the following tasks:

* removes non-canonical contigs
* adds GT and AD fields (dot or calculated based on NT, SGT, if available)
* removes tool-specific header lines

## Overview

![vmerging flowchart](docs/VARMERGE_specs.png)

## Dependencies

* [tabix 0.2.6](https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2)
* [gatk 4.2.6.1](https://gatk.broadinstitute.org)


## Usage

### Cromwell
```
java -jar cromwell.jar run variantMerging.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`reference`|String|Reference assmbly id, passed by the respective olive
`inputVcfs`|Array[Pair[File,String]]|Pairs of vcf files (SNV calls from different callers) and metadata string (producer of calls).
`tumorName`|String|Tumor id to use in vcf headers
`outputFileNamePrefix`|String|Output prefix to prefix output file names with.


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`normalName`|String?|None|Normal id to use in vcf headers, Optional


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`preprocessVcf.preprocessScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"|path to preprocessing script
`preprocessVcf.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`preprocessVcf.timeout`|Int|10|timeout in hours
`mergeVcfs.timeout`|Int|20|timeout in hours
`mergeVcfs.jobMemory`|Int|12|Allocated memory, in GB
`combineVariants.combiningScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfCombine.py"|Path to combining script
`combineVariants.jobMemory`|Int|12|memory allocated to preprocessing, in GB
`combineVariants.timeout`|Int|20|timeout in hours
`postprocessMerged.postprocessScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"|path to postprocessing script, this is the same script we use for pre-processing
`postprocessMerged.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`postprocessMerged.timeout`|Int|10|timeout in hours
`postprocessCombined.postprocessScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"|path to postprocessing script, this is the same script we use for pre-processing
`postprocessCombined.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`postprocessCombined.timeout`|Int|10|timeout in hours


### Outputs

Output | Type | Description | Labels
---|---|---|---
`mergedVcf`|File|vcf file containing all variant calls|vidarr_label: mergedVcf
`mergedIndex`|File|tabix index of the vcf file containing all variant calls|vidarr_label: mergedIndex
`combinedVcf`|File|combined vcf file containing all variant calls|vidarr_label: combinedVcf
`combinedIndex`|File|index of combined vcf file containing all variant calls|vidarr_label: combinedIndex


## Commands
 
 This section lists command(s) run by variantMerging workflow
 
### Preprocessing
 
 Detect NORMAL/TUMOR swap, impute missing fields (i.e. in case of such callers as strelka) 
 
```
  python3 PREPROCESSING_SCRIPT VCF_FILE -o VCF_FILE_BASENAME_tmp.vcf -r REFERENCE_ID
  
  bgzip -c VCF_FILE_BASENAME_tmp.vcf > VCF_FILE_BASENAME_tmp.vcf.gz
  gatk SortVcf -I VCF_FILE_BASENAME_tmp.vcf.gz 
               -R REF_FASTA 
               -O VCF_FILE_BASENAME_processed.vcf.gz
```
  
### Merging vcf files
  
This is a simple concatenation of input vcfs, there may be duplicate entries for the same call if multiple callers discover the same variant.
  
```
  gatk MergeVcfs -I INPUT_VCFS -O PREFIX_mergedVcfs.vcf.gz
  
```
  
### Combining vcf files
  
 A more complex merging: matching fields will be annotated by caller.
  
```
 
   set -euxo pipefail 
   python3 <<CODE
   import sys
   v = "~{sep=' ' inputVcfs}"
   vcfFiles = v.split()
   with open("vcf_list", 'w') as l:
       for v in vcfFiles:
           l.write(v + "\n")
   CODE
 
   python3 COMBINING_SCRIPT vcf_list -c OUTPUT_PREFIX_tmp.vcf -n ~{sep=',' inputNames}
   gatk SortVcf -I OUTPUT_PREFIX_tmp.vcf -R REFERENCE_FASTA -O OUTPUT_PREFIX_combined.vcf.gz
 
```
### Postprocessing (Name injection)
 
 The same script used for preprocessing injects names of samples into the header if argumants are passed
 
```
  set -euxo pipefail
  python3 ~{postprocessScript} ~{vcfFile} -o ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf -r ~{referenceId} -t ~{tumorName} ~{"-n " + normalName}
  bgzip -c ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf > ~{basename(vcfFile, '.vcf.gz')}.vcf.gz
  tabix -p vcf ~{basename(vcfFile, '.vcf.gz')}.vcf.gz
 
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
