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
* [bcftools 1.9](https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2)
* [gatk 4.2.6.1](https://gatk.broadinstitute.org)
* [bcbio-variation-recall 0.2.6](https://github.com/bcbio/bcbio.variation.recall/archive/refs/tags/v0.2.6.tar.gz)


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
`priorities`|Array[String]|List of workflows which produced the calls, ordered by priority
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
`resortList.timeout`|Int|20|timeout in hours
`resortList.jobMemory`|Int|12|Allocated memory, in GB
`mergeVcfsAll.timeout`|Int|20|timeout in hours
`mergeVcfsAll.jobMemory`|Int|12|Allocated memory, in GB
`combineVariantsAll.combiningScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfCombine.py"|Path to combining script
`combineVariantsAll.jobMemory`|Int|12|memory allocated to preprocessing, in GB
`combineVariantsAll.timeout`|Int|20|timeout in hours
`ensembleVariantsAll.ensembleProgram`|String|"$BCBIO_VARIATION_RECALL_ROOT/bin/bcbio-variation-recall"|Path to ensemble program
`ensembleVariantsAll.additionalParameters`|String?|None|Optional additional parameters for ensemble program
`ensembleVariantsAll.minCallers`|Int|1|variant gets recorded if minimum number of callers make the call
`ensembleVariantsAll.jobMemory`|Int|12|memory allocated to preprocessing, in GB
`ensembleVariantsAll.timeout`|Int|20|timeout in hours
`mergeVcfsPass.timeout`|Int|20|timeout in hours
`mergeVcfsPass.jobMemory`|Int|12|Allocated memory, in GB
`combineVariantsPass.combiningScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfCombine.py"|Path to combining script
`combineVariantsPass.jobMemory`|Int|12|memory allocated to preprocessing, in GB
`combineVariantsPass.timeout`|Int|20|timeout in hours
`ensembleVariantsPass.ensembleProgram`|String|"$BCBIO_VARIATION_RECALL_ROOT/bin/bcbio-variation-recall"|Path to ensemble program
`ensembleVariantsPass.additionalParameters`|String?|None|Optional additional parameters for ensemble program
`ensembleVariantsPass.minCallers`|Int|1|variant gets recorded if minimum number of callers make the call
`ensembleVariantsPass.jobMemory`|Int|12|memory allocated to preprocessing, in GB
`ensembleVariantsPass.timeout`|Int|20|timeout in hours
`postprocessMerged.postprocessScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"|path to postprocessing script, this is the same script we use for pre-processing
`postprocessMerged.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`postprocessMerged.timeout`|Int|10|timeout in hours
`postprocessCombined.postprocessScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"|path to postprocessing script, this is the same script we use for pre-processing
`postprocessCombined.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`postprocessCombined.timeout`|Int|10|timeout in hours
`postprocessEnsembled.postprocessScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"|path to postprocessing script, this is the same script we use for pre-processing
`postprocessEnsembled.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`postprocessEnsembled.timeout`|Int|10|timeout in hours
`postprocessMergedPass.postprocessScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"|path to postprocessing script, this is the same script we use for pre-processing
`postprocessMergedPass.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`postprocessMergedPass.timeout`|Int|10|timeout in hours
`postprocessCombinedPass.postprocessScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"|path to postprocessing script, this is the same script we use for pre-processing
`postprocessCombinedPass.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`postprocessCombinedPass.timeout`|Int|10|timeout in hours
`postprocessEnsembledPass.postprocessScript`|String|"$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"|path to postprocessing script, this is the same script we use for pre-processing
`postprocessEnsembledPass.jobMemory`|Int|12|memory allocated to preprocessing, in gigabytes
`postprocessEnsembledPass.timeout`|Int|10|timeout in hours


### Outputs

Output | Type | Description | Labels
---|---|---|---
`mergedVcf`|File|vcf file containing all variant calls|vidarr_label: mergedVcf
`mergedIndex`|File|tabix index of the vcf file containing all variant calls|vidarr_label: mergedIndex
`combinedVcf`|File|combined vcf file containing all variant calls|vidarr_label: combinedVcf
`combinedIndex`|File|index of combined vcf file containing all variant calls|vidarr_label: combinedIndex
`ensembledVcf`|File|endembled vcf file containing all variant calls|vidarr_label: ensembledVcf
`ensembledIndex`|File|index of ensembled vcf file containing all variant calls|vidarr_label: ensembledIndex
`mergedPassVcf`|File|vcf file containing merged PASS calls|vidarr_label: mergedPassVcf
`mergedPassIndex`|File|tabix index of the vcf file containing merged PASS calls|vidarr_label: mergedPassIndex
`combinedPassVcf`|File|combined vcf file containing combined PASS calls|vidarr_label: combinedPassVcf
`combinedPassIndex`|File|index of combined vcf file containing PASS calls|vidarr_label: combinedPassIndex
`ensembledPassVcf`|File|endembled vcf file containing PASS calls|vidarr_label: ensembledPassVcf
`ensembledPassIndex`|File|index of ensembled vcf file containing PASS calls|vidarr_label: ensembledPassIndex


## Commands
 
This section lists command(s) run by variantMerging workflow
 
* Running variantMerging
 
### Preprocessing
 
 Detect NORMAL/TUMOR swap, impute missing fields (i.e. in case of such callers as strelka) 
 A vetting script makes sure we have matching formats used across vcf, in addition making separate vcf files with only PASS calls
 
```
  set -euxo pipefail
  python3 ~{preprocessScript} ~{vcfFile} -o ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf -r ~{referenceId} 
  bgzip -c ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf > ~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz
  bcftools view -f "PASS" ~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz | bgzip -c > ~{basename(vcfFile, '.vcf.gz')}_processed_pass.vcf.gz
```
 
### reorder inputs according to priority
 
```
  python3 <<CODE
  import re
  sorted_indices = []
  unsortedNames = re.split(",",  "~{sep=',' unsortedWorkflows}")
  priorities = re.split(",", "~{sep=',' priorities}")
  unsortedFiles = re.split(",", "~{sep=',' unsortedVcfs}")
  unsortedPassFiles = re.split(",", "~{sep=',' unsortedPassVcfs}")
  sorted_indices = []
  for p in  priorities:
      if p in unsortedNames:
          print(p + "\n")
          sorted_indices.append(unsortedNames.index(p))
 
  print(sorted_indices)
  with open("~{sortedFiles}", mode='w') as out:
     out.writelines([unsortedFiles[i] + "\n" for i in sorted_indices])
  with open("~{sortedPassFiles}", mode='w') as out2:
     out2.writelines([unsortedPassFiles[j] + "\n" for j in sorted_indices])
  with open("~{sortedWorkflows}", mode='w') as out3:
     out3.writelines([unsortedNames[k] + "\n" for k in sorted_indices])
  CODE
```
 
### Merge variants with GATK (picard)
 
```
  gatk MergeVcfs -I ~{sep=" -I " inputVcfs} -O ~{outputPrefix}_mergedVcfs.vcf.gz
```
 
### Customized combining of the variants
 
 This step is custom-scripted and the produced vcf has variants annotated in a very detailed way

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

### Ensemble vcfs (combine calls using bcbio approach)
 
```
   ~{ensembleProgram} ensemble ~{outputPrefix}_ensembled.vcf.gz ~{referenceFasta} --names ~{sep=',' inputNames} --numpass ~{minCallers} ~{additionalParameters} ~{sep=' ' inputVcfs}
```
 
### Post-processing
 
```
  set -euxo pipefail
  python3 ~{postprocessScript} ~{vcfFile} -o ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf -r ~{referenceId} -t ~{tumorName} ~{"-n " + normalName}
  bgzip -c ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf > ~{basename(vcfFile, '.vcf.gz')}.vcf.gz
  tabix -p vcf ~{basename(vcfFile, '.vcf.gz')}.vcf.gz
 
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
