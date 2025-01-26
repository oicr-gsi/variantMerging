version 1.0

struct GenomeResources {
    String refModule
    String refFasta
}

workflow variantMerging {
input {
  String reference
  Array[Pair[File, String]] inputVcfs
  String tumorName
  String? normalName
  String outputFileNamePrefix
}

parameter_meta {
  reference: "Reference assmbly id, passed by the respective olive"
  inputVcfs: "Pairs of vcf files (SNV calls from different callers) and metadata string (producer of calls)."
  tumorName: "Tumor id to use in vcf headers"
  normalName: "Normal id to use in vcf headers, Optional"
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
       modules = resources[reference].refModule + " gatk/4.2.6.1 varmerge-scripts/2.1 tabix/0.2.6",
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
     modules = resources[reference].refModule + " varmerge-scripts/2.1 gatk/4.2.6.1",
     referenceFasta = resources[reference].refFasta
}

# Ensemble variant using bcbio software
call ensembleVariants {
  input:
    inputVcfs = preprocessVcf.processedVcf,
    outputPrefix = outputFileNamePrefix,
    modules = resources[reference].refModule + " bcbio-variation-recall/0.2.6",
    referenceFasta = resources[reference].refFasta
}

# Post-process vcf files, inject names into the headers for downstream compatibility
call postprocessVcf as postprocessMerged {
  input:
       vcfFile = mergeVcfs.mergedVcf,
       modules = resources[reference].refModule + " varmerge-scripts/2.1 tabix/0.2.6",
       referenceId = reference,
       tumorName = tumorName,
       normalName = normalName
}

call postprocessVcf as postprocessCombined {
  input:
       vcfFile = combineVariants.combinedVcf,
       modules = resources[reference].refModule + " varmerge-scripts/2.1 tabix/0.2.6",
       referenceId = reference,
       tumorName = tumorName,
       normalName = normalName
}

call postprocessVcf as postprocessEnsembled {
  input:
       vcfFile = ensembleVariants.ensembledVcf,
       modules = resources[reference].refModule + " varmerge-scripts/2.1 tabix/0.2.6",
       referenceId = reference,
       tumorName = tumorName,
       normalName = normalName
}

meta {
  author: "Peter Ruzanov"
  email: "pruzanov@oicr.on.ca"
  description: "VariantMerging is a workflow for combining variant calls from SNV analyses done with different callers (such as muTect2, strelka2). The workflow pre-processes input vcf files by removing non-canonical contigs, fixing fields and inferring missing values from available data. It combines calls, annotating them with caller-specific tags which allows identification of consensus variants. The workflow also uses GATK for producing merged results. In this case, all calls appear as-as. Essentially, this is a simple concatenation of the inputs.\n### Pre-processing\n\nThe script used at this step performs the following tasks:\n\n* removes non-canonical contigs\n* adds GT and AD fields (dot or calculated based on NT, SGT, if available)\n* removes tool-specific header lines\n\n## Overview\n\n![vmerging flowchart](docs/VARMERGE_specs.png)"
  dependencies: [
      {
        name: "tabix/0.2.6",
        url: "https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2"
      },
      {
        name: "gatk/4.2.6.1",
        url: "https://gatk.broadinstitute.org"
      },
      {
        name: "bcbio-variation-recall/0.2.6",
        url: "https://github.com/bcbio/bcbio.variation.recall/archive/refs/tags/v0.2.6.tar.gz"
      }
    ]
    output_meta: {
    mergedVcf: {
        description: "vcf file containing all variant calls",
        vidarr_label: "mergedVcf"
    },
    mergedIndex: {
        description: "tabix index of the vcf file containing all variant calls",
        vidarr_label: "mergedIndex"
    },
    combinedVcf: {
        description: "combined vcf file containing all variant calls",
        vidarr_label: "combinedVcf"
    },
    combinedIndex: {
        description: "index of combined vcf file containing all variant calls",
        vidarr_label: "combinedIndex"
    },
    ensembledVcf: {
        description: "endembled vcf file containing all variant calls",
        vidarr_label: "ensembledVcf"
    },
    ensembledIndex: {
        description: "index of ensembled vcf file containing all variant calls",
        vidarr_label: "ensembledIndex"
    }
}
}

output {
  File mergedVcf = postprocessMerged.processedVcf
  File mergedIndex = postprocessMerged.processedIndex
  File combinedVcf = postprocessCombined.processedVcf
  File combinedIndex = postprocessCombined.processedIndex
  File ensembledVcf = postprocessEnsembled.processedVcf
  File ensembledIndex = postprocessEnsembled.processedIndex
}
}

# =================================
#  1 of 5: preprocess a vcf file 
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
 bgzip -c ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf > ~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz
 tabix -p vcf ~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz
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
#  2 of 5: concat files, all variants concatenated,
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
#  3 of 5: merge files with CombineVariants, merging matching fields
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

# ==================================================================
#  4 of 5: merge files with bcbio tool bcbio-variation-recall
#  with priority defined by the order in the input array
# ==================================================================
task ensembleVariants {
input {
 Array[File] inputVcfs
 String outputPrefix
 String modules 
 String ensembleProgram = "$BCBIO_VARIATION_RECALL_ROOT/bin/bcbio-variation-recall"
 String referenceFasta
 String? additionalParameters
 Int jobMemory = 12
 Int timeout = 20
}

parameter_meta {
 inputVcfs: "array of input vcf files"
 outputPrefix: "prefix for output file"
 modules: "modules for running preprocessing"
 ensembleProgram: "Path to ensemble program"
 referenceFasta: "path to the reference FASTA file"
 additionalParameters: "Optional additional parameters for ensemble program"
 jobMemory: "memory allocated to preprocessing, in GB"
 timeout: "timeout in hours"
}

command <<<
  ~{ensembleProgram} ensemble ~{outputPrefix}_ensembled.vcf.gz ~{referenceFasta} ~{additionalParameters} ~{sep=' ' inputVcfs}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File ensembledVcf = "~{outputPrefix}_ensembled.vcf.gz"
  File ensembledIndex = "~{outputPrefix}_ensembled.vcf.gz.tbi"
}
}

# =======================================================================
#  5 of 5: post-preprocess a vcf file (this is mainly for name injection)
# =======================================================================
task postprocessVcf {
input {
 File vcfFile
 String referenceId
 String postprocessScript = "$VARMERGE_SCRIPTS_ROOT/bin/vcfVetting.py"
 String modules
 String tumorName
 String? normalName
 Int jobMemory = 12
 Int timeout = 10
}

parameter_meta {
 vcfFile: "path to the input vcf file"
 referenceId: "String that shows the id of the reference assembly"
 postprocessScript: "path to postprocessing script, this is the same script we use for pre-processing"
 tumorName: "Tumor name to use in vcf header"
 normalName: "Normal name to use in vcf header, Optional"
 modules: "modules for running preprocessing"
 jobMemory: "memory allocated to preprocessing, in gigabytes"
 timeout: "timeout in hours"
}

command <<<
 set -euxo pipefail
 python3 ~{postprocessScript} ~{vcfFile} -o ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf -r ~{referenceId} -t ~{tumorName} ~{"-n " + normalName}
 bgzip -c ~{basename(vcfFile, '.vcf.gz')}_tmp.vcf > ~{basename(vcfFile, '.vcf.gz')}.vcf.gz
 tabix -p vcf ~{basename(vcfFile, '.vcf.gz')}.vcf.gz
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File processedVcf   = "~{basename(vcfFile, '.vcf.gz')}.vcf.gz"
  File processedIndex = "~{basename(vcfFile, '.vcf.gz')}.vcf.gz.tbi"
}
}
