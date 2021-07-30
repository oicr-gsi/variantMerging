# variantMerging

VariantMerging 2.0, a workflow for combining variant calls from SNV analyses done with different callers

### Pre-processing

The script used at this step performs the following tasks:

* removes non-canonical contigs
* adds GT and AD fields (dot or calculated based on NT, SGT, if available)
* removes tool-specific header lines

## Overview

![vmerging flowchart](docs/VARMERGE_specs.png)

## Dependencies

* [java 8](https://github.com/AdoptOpenJDK/openjdk8-upstream-binaries/releases/download/jdk8u222-b10/OpenJDK8U-jdk_x64_linux_8u222b10.tar.gz)
* [tabix 0.2.6](https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2)
* [gatk 4.1.7.0, gatk 3.6.0](https://gatk.broadinstitute.org)


## Usage

### Cromwell
```
java -jar cromwell.jar run variantMerging.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputVcfs`|Array[Pair[File,String]]|Pairs of vcf files (SNV calls from different callers) and metadata string (producer of calls).
`preprocessVcf.referenceId`|String|String that shows the id of the reference assembly
`preprocessVcf.referenceFasta`|String|path to the reference FASTA file
`preprocessVcf.preprocessScript`|String|path to preprocessing script
`preprocessVcf.modules`|String|modules for running preprocessing
`mergeVcfs.modules`|String|modules for this task
`combineVariants.referenceFasta`|String|path to the reference FASTA file
`combineVariants.modules`|String|modules for running preprocessing
`combineVariants.priority`|String|Comma-separated list defining priority of workflows when combining variants


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|""|Output prefix to prefix output file names with.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`preprocessVcf.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`preprocessVcf.timeout`|Int|10|timeout in hours
`mergeVcfs.timeout`|Int|20|timeout in hours
`mergeVcfs.jobMemory`|Int|12|Allocated memory, in GB
`combineVariants.jobMemory`|Int|12|memory allocated to preprocessing, in GB
`combineVariants.timeout`|Int|20|timeout in hours


### Outputs

Output | Type | Description
---|---|---
`mergedVcf`|File|vcf file containing all structural variant calls
`mergedIndex`|File|tabix index of the vcf file containing all structural variant calls
`combinedVcf`|File|filtered vcf file containing all structural variant calls
`combinedIndex`|File|tabix index of the filtered vcf file containing all structural variant calls


## Commands
 
 This section lists commands run by variantMerging workflow
 
 ### Preprocessing
 
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
 
 A more complex merging with GATK CombineVariants: depending on priority assigned to the callers matching fields will be ranked according this settings and only the values from the caller with highest priority will be used.
 
 
 ```
   ...
   
   Embedded Python code runs the CombineVariants command:
 
   java -Xmx[JOB_MEMORY]G -jar GenomeAnalysisTK.jar
        -T CombineVariants INPUTS
        -R EF_FASTA
        -o PREFIX_combined.vcf.gz
        -genotypeMergeOptions PRIORITIZE
        -priority PRIORITY
 
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
