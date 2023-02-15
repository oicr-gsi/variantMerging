version 1.0

workflow variantMerging {
input {
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
# Merge using MergeVcfs (Picard) 
call mergeVcfs { input: inputVcfs = preprocessVcf.processedVcf, outputPrefix = sampleID }

# Combine using DISCVR-Seq Toolkit
call combineVariants {input: inputVcfs = preprocessVcf.processedVcf, inputIndexes = preprocessVcf.processedIndex, workflows = preprocessVcf.prodWorkflow, outputPrefix = sampleID }

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "VariantMerging 2.1, a workflow for combining variant calls from SNV analyses done with different callers\n### Pre-processing\n\nThe script used at this step performs the following tasks:\n\n* removes non-canonical contigs\n* adds GT and AD fields (dot or calculated based on NT, SGT, if available)\n* removes tool-specific header lines\n\n## Overview\n\n![vmerging flowchart](docs/VARMERGE_specs.png)"
  dependencies: [
     {
        name: "java/9",
        url: "https://github.com/AdoptOpenJDK/openjdk9-binaries/releases/download/jdk-9%2B181/OpenJDK9U-jdk_x64_linux_hotspot_9_181.tar.gz"
      },
      {
        name: "tabix/0.2.6",
        url: "https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2"
      },
      {
        name: "gatk/4.2.6.1",
        url: "https://gatk.broadinstitute.org"
      },
      {
        name: "DISCVR-Seq Toolkit/1.3.21",
        url: "https://bimberlab.github.io/DISCVRSeq"
      }
    ]
    output_meta: {
      mergedVcf: "vcf file containing all variant calls",
      mergedIndex: "tabix index of the vcf file containing all variant calls",
      combinedVcf: "combined vcf file containing all variant calls",
      combinedIndex: "index of combined vcf file containing all variant calls"
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
 String modules = "gatk/4.2.6.1 varmerge-scripts/1.6 tabix/0.2.6"
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
 String modules = "gatk/4.2.6.1 tabix/0.2.6"
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
 String modules = "discvrseq/1.3.21"
 String priority
 String dscrvToolsJar = "$DISCVRSEQ_ROOT/bin/DISCVRSeq-1.3.21.jar"
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
 dscrvToolsJar: "DISCVR tools JAR"
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

  dscvrCommand  = "$JAVA_ROOT/bin/java -Xmx" + str(javaMemory) + "G -jar ~{dscrvToolsJar} "
  dscvrCommand += "MergeVcfsAndGenotypes "
  dscvrCommand += " ".join(inputStrings)
  dscvrCommand += " -R ~{referenceFasta} "
  dscvrCommand += "-O ~{outputPrefix}_combined.vcf.gz "
  dscvrCommand += "--genotypeMergeOption PRIORITIZE "
  dscvrCommand += "-priority " + priority
  dscvrCommand += " 2>&1"

  result_output = subprocess.run(dscvrCommand, shell=True)
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
