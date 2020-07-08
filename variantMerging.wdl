version 1.0

workflow variantMerging {
input {
    # Normally we need only an array of vcf.gz files, first one gets the highest
    # priority if there are conflicting values for matching fields
    Array[Pair[File, String]] inputVcfs
    String? outputFileNamePrefix = ""
}

String? sampleID = if outputFileNamePrefix=="" then basename(inputVcfs[0].left), ".vcf.gz") else outputFileNamePrefix

parameter_meta {
  inputVcfs: "Pairs of vcf files (SNV calls from different callers) and metadata string (producer of calls)."
  outputFileNamePrefix: "Output prefix to prefix output file names with."
}

# Preprocess vcf files
scatter (v in inputVcfs) {
  call preprocessVcf { input: vcfFile = v }
}

# Do merging (two ways)

# Post-process not really needed?
meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "VariantMerging 2.0"
  
}

output {
  File mergedVcf   = mergeVcfs.infoFile
  File combinedVcf = combineVariants.regionFile
}

}

# ==============================================
#  1 of 3: preprocess a vcf file (and validate)
# ==============================================
task preprocessVcf {
input {
 File vcfFile
 String modules = ""
 Int jobMemory  = 12
 Int timeout    = 10
}

parameter_meta {
 inputTumor: "Input .bam file for analysis sample"
 inputNormal: "Optional input .bam file for control sample"
}

command <<<
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File processedVcf = "basename(inputVcfs[0], ".vcf.gz")_processed.vcf.gz"
}
}

# ==================================================
#  2a of 3: concat files, all variants concatenated,
#  same variant may appear multiple times
# ==================================================
task mergeVcfs {
input {
 Array[File] inputVcfs
 Int timeout
 Int jobMemory
 String modules = "gatk/4.1.7.0"
}

parameter_meta {
 inputVcfs: "Array of vcf files to merge"
 inputTumor: "Input .bam file for analysis sample"
 inputNormal: "Optional input .bam file for control sample"
}

command <<<
 gatk MergeVcfs ~{sep="-I ", inputVcfs} -O ~{outputPrefix}_mergedVcfs.vcf.gz
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File mergedVcf = "~{outputPrefix}_mergedVcfs.vcf.gz"
}
}

# ==================================================================
#  2b of 3: merge files with CombineVariants, merging matching fields
#  with priority defined by the order in the input array
# ==================================================================
task combineVariants {
input {
 Array[File] inputs
 Array[String] workflows
 String reference_fasta
 String outputPrefix
 String modules
 Int jobMemory
 Int timeout
}


parameter_meta {
 inputTumor: "Input .bam file for analysis sample"
 inputNormal: "Optional input .bam file for control sample"
}

# perhaps python embedded code is better here
command <<<
  $JAVA_HOME/bin/java -jar $GATK_ROOT/GenomeAnalysisTK.jar \
                      -T CombineVariants \
                      ~{sep=" " inputStrings} \ # strings "--variant:workflow FilePath" ???
                      -R ~{reference_fasta} \
                      -o ~{outputPrefix}_combined.vcf.gz \
                      -genotypeMergeOptions PRIORITIZE \
                      -priority ~{sep="," workflows}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File combinedVcf = "~{outputPrefix}_combined.vcf.gz"
}
}

# ===================================================
#  3 of 3: post-process files, merge columns
# ===================================================
task postprocessVcf {
input {

}

parameter_meta {
 inputTumor: "Input .bam file for analysis sample"
 inputNormal: "Optional input .bam file for control sample"
}

command <<<
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File? ratioBedGraph = "~{sampleID}_processed.vcf.gz"
}
}
