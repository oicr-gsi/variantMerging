version 1.0

workflow variantMerging {
input {
    # Normally we need only an array of vcf.gz files, first one gets the highest
    # priority if there are conflicting values for matching fields
    Array[Pair[File, String]] inputVcfs
    String outputFileNamePrefix = ""
}

String sampleID = if outputFileNamePrefix=="" then basename(inputVcfs[0].left, ".vcf.gz") else outputFileNamePrefix

parameter_meta {
  inputVcfs: "Pairs of vcf files (SNV calls from different callers) and metadata string (producer of calls)."
  outputFileNamePrefix: "Output prefix to prefix output file names with."
}

# Preprocess vcf files
scatter (v in inputVcfs) {
  call preprocessVcf { input: vcfFile = v.left, producerWorkflow = v.right }
}

# Do merging (two ways)
call mergeVcfs { input: inputVcfs = preprocessVcf.processedVcf, outputPrefix = sampleID }

call combineVariants {input: inputVcfs = preprocessVcf.processedVcf, inputIndexes = preprocessVcf.processedIndex, workflows = preprocessVcf.prodWorkflow, outputPrefix = sampleID }

# Post-process not really needed?
meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "VariantMerging 2.0, a workflow for combining variant calls from SNV analyses done with different callers"
  dependencies: [
     {
        name: "java/8",
        url: "https://github.com/AdoptOpenJDK/openjdk8-upstream-binaries/releases/download/jdk8u222-b10/OpenJDK8U-jdk_x64_linux_8u222b10.tar.gz"
      },
      {
        name: "tabix/0.2.6",
        url: "https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2"
      },
      {
        name: "gatk/4.1.7.0, gatk/3.6.0",
        url: "https://gatk.broadinstitute.org"
      }
    ]
    output_meta: {
      mergedVcf: "vcf file containing all structural variant calls",
      mergedIndex: "tabix index of the vcf file containing all structural variant calls",
      combinedVcf: "filtered vcf file containing all structural variant calls",
      combinedIndex: "tabix index of the filtered vcf file containing all structural variant calls"
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
 String preprocessScript
 String modules
 Int jobMemory  = 12
 Int timeout    = 10
}

parameter_meta {
 vcfFile: "path to the input vcf file"
 producerWorkflow: "workflow name that produced the vcf"
 referenceId: "String that shows the id of the reference assembly"
 preprocessScript: "path to preprocessing script"
 modules: "modules for running preprocessing"
 jobMemory: "memory allocated to preprocessing, in gigabytes"
 timeout: "timeout in hours"
}

command <<<
 set -euxo pipefail
 python3 ~{preprocessScript} ~{vcfFile} -o ~{basename(vcfFile, '.vcf.gz')}_processed.vcf -r ~{referenceId}
 bgzip -c ~{basename(vcfFile, '.vcf.gz')}_processed.vcf > ~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz
 tabix -p vcf ~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File processedVcf = "~{basename(vcfFile, '.vcf.gz')}_processed.vcf.gz"
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
 Array[File] inputIndexes
 Array[String] workflows
 String referenceFasta
 String outputPrefix
 String modules
 String priority
 Int jobMemory = 12
 Int timeout = 20
}

parameter_meta {
 inputVcfs: "array of input vcf files"
 inputIndexes: "array of tabix indexes for vcf files"
 workflows: "array of ids of producer workflows"
 referenceFasta: "path to the reference FASTA file"
 outputPrefix: "prefix for output file"
 modules: "modules for running preprocessing"
 priority: "Comma-separated list defining priority of workflows when combining variants"
 jobMemory: "memory allocated to preprocessing, in GB"
 timeout: "timeout in hours"
}

command <<<
  python3<<CODE
  import subprocess
  import sys
  inputStrings = []
  v = "~{sep=' ' inputVcfs}"
  vcfFiles = v.split()
  w = "~{sep=' ' workflows}"
  workflowIds = w.split()
  priority = "~{priority}"
  
  if len(vcfFiles) != len(workflowIds):
      print("The arrays with input files and their respective workflow names are not of equal size!")
  else:
      for f in range(0, len(vcfFiles)):
          inputStrings.append("--variant:" + workflowIds[f] + " " + vcfFiles[f])

  javaMemory = ~{jobMemory} - 6 
  gatkCommand  = "$JAVA_ROOT/bin/java -Xmx" + str(javaMemory) + "G -jar $GATK_ROOT/GenomeAnalysisTK.jar "
  gatkCommand += "-T CombineVariants "
  gatkCommand += " ".join(inputStrings)
  gatkCommand += " -R ~{referenceFasta} "
  gatkCommand += "-o ~{outputPrefix}_combined.vcf.gz "
  gatkCommand += "-genotypeMergeOptions PRIORITIZE "
  gatkCommand += "-priority " + priority
  gatkCommand += " 2>&1"

  result_output = subprocess.run(gatkCommand, shell=True)
  sys.exit(result_output.returncode)
  CODE
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
