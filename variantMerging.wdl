version 1.0

struct GenomeResources {
    String refModule
    String refFasta
}

workflow variantMerging {
input {
  String reference
  Array[Pair[File, String]] inputVcfs
  String outputFileNamePrefix
}

parameter_meta {
  reference: "Reference assmbly id, passed by the respective olive"
  inputVcfs: "Pairs of vcf files (SNV calls from different callers) and metadata string (producer of calls)."
  outputFileNamePrefix: "Output prefix to prefix output file names with."
}

Map[String,GenomeResources] resources = {
  "hg19": {
    "refModule": "hg19/p13", 
    "refFasta": "$HG19_ROOT/hg19_random.fa"
  },
  "hg38": {
    "refModule": "hg38/p12",
    "refFasta": "$HG38_ROOT/hg38_random.fa"
  },
  "mm10": {
    "refModule": "mm10/p6",
    "refFasta": "$MM10_ROOT/mm10.fa"
  }
}

# Preprocess vcf files
scatter (v in inputVcfs) {
  call preprocessVcf {
    input: 
       vcfFile = v.left, 
       producerWorkflow = v.right, 
       modules = resources[reference].refModule + " gatk/4.2.6.1 varmerge-scripts/2.0 tabix/0.2.6", 
       referenceId = reference,
       referenceFasta = resources[reference].refFasta
  }
}

# Do merging (two ways)
# Merge using MergeVcfs (Picard) 
call mergeVcfs {
  input:
     inputVcfs = preprocessVcf.processedVcf,
     outputPrefix = outputFileNamePrefix,
     modules = "gatk/4.2.6.1 tabix/0.2.6"
}

# Combine using Custom script
call combineVariants {
  input: 
     inputVcfs = preprocessVcf.processedVcf,
     inputNames = preprocessVcf.prodWorkflow,
     outputPrefix = outputFileNamePrefix,
     modules = resources[reference].refModule + " varmerge-scripts/2.0 gatk/4.2.6.1",
     referenceFasta = resources[reference].refFasta
}

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "VariantMerging 2.1, a workflow for combining variant calls from SNV analyses done with different callers\n### Pre-processing\n\nThe script used at this step performs the following tasks:\n\n* removes non-canonical contigs\n* adds GT and AD fields (dot or calculated based on NT, SGT, if available)\n* removes tool-specific header lines\n\n## Overview\n\n![vmerging flowchart](docs/VARMERGE_specs.png)"
  dependencies: [
      {
        name: "tabix/0.2.6",
        url: "https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2"
      },
      {
        name: "gatk/4.2.6.1",
        url: "https://gatk.broadinstitute.org"
      }
    ]
    output_meta: {
      mergedVcf: "vcf file containing all variant calls",
      mergedIndex: "tabix index of the vcf file containing all variant calls",
      combinedVcf: "combined vcf file containing all variant calls",
      combinedIndex: "index of combined vcf file containing all variant calls",
      postprocessedVcf: "post-processed combined vcf file with updated set field for overlapping calls",
      postprocessedIndex: "index of post-processed vcf file"
    }
}

output {
  File mergedVcf = mergeVcfs.mergedVcf
  File mergedIndex = mergeVcfs.mergedIndex
  File combinedVcf = combineVariants.combinedVcf
  File combinedIndex = combineVariants.combinedIndex
}

}

# =================================
#  1 of 3: preprocess a vcf file 
# =================================
task preprocessVcf {
input {
 File vcfFile
 String producerWorkflow
 String referenceId
 String referenceFasta
 String preprocessScript = "$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"
 String modules
 Int jobMemory = 12
 Int timeout = 10
}

parameter_meta {
 vcfFile: "path to the input vcf file"
 producerWorkflow: "workflow name that produced the vcf"
 referenceId: "String that shows the id of the reference assembly"
 referenceFasta: "path to the reference FASTA file"
 preprocessScript: "path to preprocessing script"
 modules: "modules for running preprocessing"
 jobMemory: "memory allocated to preprocessing, in gigabytes"
 timeout: "timeout in hours"
}

command <<<
 set -euxo pipefail
 python3 ~{preprocessScript} ~{vcfFile} -o ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf -r ~{referenceId}
 bgzip -c ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf > ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf.gz
 gatk SortVcf -I ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf.gz -R ~{referenceFasta} -O ~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File processedVcf   = "~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz"
  File processedIndex = "~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz.tbi"
  String prodWorkflow = producerWorkflow
}
}

# ==================================================
#  2 of 3: concat files, all variants concatenated,
#  same variant may appear multiple times
# ==================================================
task mergeVcfs {
input {
 Array[File] inputVcfs
 String outputPrefix
 Int timeout = 20
 Int jobMemory = 12
 String modules
}

parameter_meta {
 inputVcfs: "Array of vcf files to merge"
 outputPrefix: "prefix for output file"
 timeout: "timeout in hours" 
 jobMemory: "Allocated memory, in GB"
 modules: "modules for this task"
}

command <<<
 gatk MergeVcfs -I ~{sep=" -I " inputVcfs} -O ~{outputPrefix}_mergedVcfs.vcf.gz
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File mergedVcf = "~{outputPrefix}_mergedVcfs.vcf.gz"
  File mergedIndex = "~{outputPrefix}_mergedVcfs.vcf.gz.tbi"
}
}

# ==================================================================
#  3 of 3: merge files with CombineVariants, merging matching fields
#  with priority defined by the order in the input array
# ==================================================================
task combineVariants {
input {
 Array[File] inputVcfs
 Array[String] inputNames
 String outputPrefix
 String modules
 String combiningScript = "$VARMERGE_SCRIPTS_ROOT/bin/vcfCombine.py"
 String referenceFasta
 Int jobMemory = 12
 Int timeout = 20
}

parameter_meta {
 inputVcfs: "array of input vcf files"
 inputNames: "Array of names, in the same order as vcf files"
 outputPrefix: "prefix for output file"
 modules: "modules for running preprocessing"
 combiningScript: "Path to combining script"
 referenceFasta: "path to the reference FASTA file"
 jobMemory: "memory allocated to preprocessing, in GB"
 timeout: "timeout in hours"
}

command <<<
  set -euxo pipefail 
  python3 <<CODE
  import sys
  v = "~{sep=' ' inputVcfs}"
  vcfFiles = v.split()
  with open("vcf_list", 'w') as l:
      for v in vcfFiles:
          l.write(v + "\n")
  CODE

  python3 ~{combiningScript} vcf_list -c ~{outputPrefix}_tmp.vcf -n ~{sep=',' inputNames}
  gatk SortVcf -I ~{outputPrefix}_tmp.vcf -R ~{referenceFasta} -O ~{outputPrefix}_combined.vcf.gz
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File combinedVcf = "~{outputPrefix}_combined.vcf.gz"
  File combinedIndex = "~{outputPrefix}_combined.vcf.gz.tbi"
}
}

