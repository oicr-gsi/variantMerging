package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.SemanticWorkflow;
import ca.on.oicr.pde.utilities.workflows.jobfactory.PicardTools;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import net.sourceforge.seqware.pipeline.workflowV2.model.Workflow;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * Basic steps: Vet the vcfs Sort with Picard Combine with GATK Select intersect
 * Post-process
 *
 * Inputs - there should be only two files, the order defines priority
 *
 */
public class VariantMergingWorkflow extends SemanticWorkflow {

    private String dataDir, tmpDir;
    private String java;

    private String picard_dir;
    private String gatk_dir;
    private String identifier;
    private String alignerName;
    private Workflow wf;
    private String tabixVersion;
    private final static String VCF_GZIPPED_METATYPE = "application/vcf-4-gzip";
    private final static String TBI_INDEX_METATYPE   = "application/tbi";
    private boolean manualOutput = false;
    private SqwFile[] inputVcfFiles;
    private ArrayList<String> inputBasenames;
    private ArrayList<String> sources;

    //References
    private String genomeFile;
    private String genomeDictionary;
    private String queue;
    private String binDir;

    private String bundledJRE;
    //Ontology-related variables
    private static final String EDAM = "EDAM";
    private static final Map<String, Set<String>> cvTerms;

    //TODO modify this
    static {
        cvTerms = new HashMap<String, Set<String>>();
        cvTerms.put(EDAM, new HashSet<String>(
                          Arrays.asList("VCF", "SNP", "tabix",
                                        "Variant Calling", "Sequence variation analysis")));
    }

    /**
     * Here we need to register all of our CV terms for attaching to result
     * files (Formats, Data, Processes)
     *
     * @return myTerms
     */
    @Override
    protected Map<String, Set<String>> getTerms() {
        Map<String, Set<String>> myTerms = new HashMap<String, Set<String>>();
        myTerms.putAll(cvTerms);
        return myTerms;
    }

    /**
     * Launched from setupDirectory();
     */
    private void VariantMergingWorkflow() {
        identifier = getProperty("identifier");
        java = getProperty("java");
        dataDir = "data";
        tmpDir = getProperty("tmp_dir");
        manualOutput = Boolean.parseBoolean(getProperty("manual_output"));

        //Picard
        picard_dir = getProperty("picard_dir");
        wf = this.getWorkflow();
        binDir = this.getWorkflowBaseDir() + "/bin/";

        alignerName = getOptionalProperty("aligner_name", "");
    }

    @Override
    public Map<String, SqwFile> setupFiles() {

        //References
        this.genomeFile = getProperty("ref_fasta");
        String derivedDictionary = this.genomeFile.substring(0, this.genomeFile.lastIndexOf(".fa")) + ".dict";
        this.genomeDictionary = getOptionalProperty("ref_dictionary", derivedDictionary);

        //GATK
        this.queue = getOptionalProperty("queue", "");
        this.gatk_dir = getOptionalProperty("gatk_dir", binDir);

        if (getProperty("tabix_version") == null) {
            Logger.getLogger(VariantMergingWorkflow.class.getName()).log(Level.SEVERE, "tabix_version is not set, we need it to call tabix correctly");
            return (null);
        } else {
            this.tabixVersion = getProperty("tabix_version");
        }

        if (getProperty("gatk_java") == null) {
            Logger.getLogger(VariantMergingWorkflow.class.getName()).log(Level.SEVERE, "gatk_java is not set, we need it to run GATK since it requires 1.7 java and we cannot rely on the defaul");
            return (null);
        } else {
            this.java = getProperty("gatk_java");
        }
        List<String> inputFilesList = Arrays.asList(StringUtils.split(getProperty("input_files"), ","));
        Set<String> inputFilesSet = new HashSet<>(inputFilesList);

        List<String> inputSourcesList = Arrays.asList(StringUtils.split(getProperty("input_sources"), ","));
        Set<String> inputSourcesSet = new HashSet<>(inputFilesList);

        if (inputFilesList.size() != inputFilesSet.size() || inputSourcesList.size() != inputSourcesSet.size()) {
            throw new RuntimeException("Duplicate sources detected in input_sources");
        }

        if (inputFilesList.size() != 2) {
            throw new RuntimeException("There must be two input files");
        }

        this.sources = new ArrayList<String>(inputSourcesList);

        return this.getFiles();
    }

    @Override
    public void setupDirectory() {
        try {
            VariantMergingWorkflow();
            this.addDirectory(dataDir);
            this.addDirectory(tmpDir);
            if (!dataDir.endsWith("/")) {
                dataDir += "/";
            }
            if (!tmpDir.endsWith("/")) {
                tmpDir += "/";
            }

            this.inputVcfFiles = provisionInputFiles("input_files");
            this.inputBasenames = new ArrayList<String>(this.inputVcfFiles.length);
            for (int i = 0; i < this.inputVcfFiles.length; i++) {
                this.inputBasenames.add(i, this.makeBasename(this.inputBasenames.get(i)));
            }
        } catch (Exception e) {
            e.printStackTrace();
            Log.error("Error in setupDirectory", e);
            ret.setReturnValue(ReturnValue.INVALIDFILE);
        }
    }

    @Override
    public void buildWorkflow() {

        // use these two strings to construct the file names for the intermediate
        // files. This way the file name describes the operations performed on it
        String operationsOnFile = "";
        String outputFile;
        String inputFile;

        Map<String, String> inputsToCombine = new HashMap<String, String>();
        List<Job> upstreamJobs = new ArrayList<Job>();

        /**
         * ==== Vet AND Sort vcf files here ====
         */
        for (int i = 0; i < this.inputBasenames.size(); i++) {
            // vet the vcfs
            operationsOnFile = ".vetted.";

            inputFile = this.inputVcfFiles[i].getProvisionedPath();
            outputFile = this.inputBasenames.get(i) + operationsOnFile;

            Job vetVcf = this.vetVcf(inputFile, this.dataDir + outputFile + ".vcf");
            vetVcf.addFile(this.inputVcfFiles[i]);

            // update seq dictionary / sort
            operationsOnFile += "sorted.";

            inputFile = outputFile;
            outputFile = this.inputBasenames.get(i) + operationsOnFile;

            Job sortVcf = this.sortVcf(this.dataDir + inputFile + ".vcf", this.dataDir + outputFile + ".vcf");
            inputsToCombine.put(this.sources.get(i), this.dataDir + outputFile + ".vcf");
            sortVcf.addParent(vetVcf);
            upstreamJobs.add(sortVcf);
        }

        /**
         * ==== Combine and Subset (extract intersect) here
         */
        for (int p = 0; p < this.sources.size(); p++) {

            // Combine based on priority
            String operationsOnMergedFile = operationsOnFile + "combined.";
            String priorityOrder = p == 0 ? this.sources.get(0) + "." + this.sources.get(1)
                                          : this.sources.get(1) + "." + this.sources.get(0);
            outputFile = identifier + operationsOnMergedFile + priorityOrder;
            Job combineVcfs = this.combineVcf(inputsToCombine, this.dataDir + outputFile + ".vcf", 0);

            for (Job j : upstreamJobs) {
                combineVcfs.addParent(j);
            }

            // produce intersect (for both overlaps)
            operationsOnMergedFile += "intersect.";
            inputFile = outputFile;
            outputFile = identifier + operationsOnMergedFile + priorityOrder;
            Job subsetVcf = this.subsetVcf(this.dataDir + inputFile + ".vcf", this.dataDir + outputFile + ".vcf");

            subsetVcf.addParent(combineVcfs);
            
            //post-process vcf (bgzip and index)
            Job postprocessVcf = this.prepareVcf(this.dataDir + outputFile + ".vcf");
            postprocessVcf.addParent(subsetVcf);
            
            // Annotate and provision final vcf and its index
            SqwFile finalVcf    = this.createOutputFile(this.dataDir + outputFile + ".vcf.gz", VCF_GZIPPED_METATYPE, manualOutput);
            SqwFile finalVcfTbi = this.createOutputFile(this.dataDir + outputFile + ".vcf.gz.tbi", TBI_INDEX_METATYPE, manualOutput);
            
            if (!this.alignerName.isEmpty()) {
                finalVcf.getAnnotations().put("aligner", this.alignerName);
                finalVcfTbi.getAnnotations().put("aligner", this.alignerName);
            }
            
            this.attachCVterms(finalVcf, EDAM, "VCF,SNP,Variant Calling,Sequence variation analysis");
            this.attachCVterms(finalVcfTbi, EDAM, "tabix,Sequence variation analysis");
            
            postprocessVcf.addFile(finalVcf);                     
            postprocessVcf.addFile(finalVcfTbi);
            
        }
    }

    //=======================Jobs as functions===================
    /**
     *
     * @param inputFile
     * @param outputFile
     * @return
     */
    protected Job vetVcf(String inputFile, String outputFile) {
        Job jobVcfPrep = wf.createBashJob("preprocess_vcf");
        jobVcfPrep.setCommand(getWorkflowBaseDir() + "/dependencies/vcf_vetting.pl "
                + inputFile + " > "
                + outputFile);

        jobVcfPrep.setMaxMemory("4000");

        if (!this.queue.isEmpty()) {
            jobVcfPrep.setQueue(this.queue);
        }

        return jobVcfPrep;
    }

    /**
     *
     * @param inputFile
     * @param outputFile
     * @return
     */
    protected Job sortVcf(String inputFile, String outputFile) {

        Job jobUpdate = wf.createBashJob("update_dictionary");
        String tempFile = inputFile + ".updated";

        jobUpdate.setCommand(getWorkflowBaseDir() + "/bin/" + this.bundledJRE + "/bin/java"
                + " -Xmx" + this.getProperty("picard_memory") + "M"
                + " -jar " + this.picard_dir + "/picard.jar UpdateVcfSequenceDictionary "
                + " I=" + inputFile
                + " O=" + tempFile
                + " SEQUENCE_DICTIONARY=" + this.genomeDictionary
                + " CREATE_INDEX=false");

        Integer mbRam = Integer.valueOf(this.getProperty("picard_memory"));
        jobUpdate.setMaxMemory(Integer.toString(mbRam + 2000));

        Job jobSort = wf.createBashJob("sort_variants");

        jobSort.setCommand(getWorkflowBaseDir() + "/bin/" + this.bundledJRE + "/bin/java"
                + " -Xmx" + this.getProperty("picard_memory") + "M"
                + " -jar " + this.picard_dir + "/picard.jar SortVcf"
                + " -R " + this.genomeFile
                + " I=" + tempFile
                + " O=" + outputFile
                + " CREATE_INDEX=false");

        jobSort.addParent(jobUpdate);
        jobSort.setMaxMemory(Integer.toString(mbRam + 2000));

        if (!this.queue.isEmpty()) {
            jobSort.setQueue(this.queue);
        }

        return jobSort;
    }

    /**
     * combine the vcf files into one according to priorities
     *
     * @param inputs
     * @param outputFile
     * @param priorityIndex
     * @return
     */
    protected Job combineVcf(Map<String, String> inputs, String outputFile, int priorityIndex) {
        Job jobGATK = this.wf.createBashJob("combine_calls");
        StringBuilder gatkCommand = new StringBuilder();
        int otherIndex = priorityIndex == 0 ? 1 : 0;

        gatkCommand.append(java)
                .append(" -Xmx4g -Djava.io.tmpdir=").append(this.tmpDir)
                .append(" -jar ").append(this.gatk_dir).append("/GenomeAnalysisTK.jar")
                .append(" -T CombineVariants ")
                .append(" --variant:").append(this.sources.get(priorityIndex)).append(" ").append(inputs.get(this.sources.get(priorityIndex)))
                .append(" --variant:").append(this.sources.get(otherIndex)).append(" ").append(inputs.get(this.sources.get(otherIndex)))
                .append(" -R ").append(this.genomeFile)
                .append(" -o ").append(outputFile)
                .append(" -genotypeMergeOptions PRIORITIZE")
                .append(" -priority ").append(this.sources.get(priorityIndex)).append(",").append(this.sources.get(otherIndex));

        jobGATK.setCommand(gatkCommand.toString());
        jobGATK.setMaxMemory("6000");

        if (!this.queue.isEmpty()) {
            jobGATK.setQueue(this.queue);
        }

        return jobGATK;
    }

    /**
     * Produce intersect subset of the calls
     *
     * @param inputFile 
     * @param outputFile
     * @return
     */
    protected Job subsetVcf(String inputFile, String outputFile) {
        Job jobGATK = this.wf.createBashJob("subset_calls");
        StringBuilder gatkCommand = new StringBuilder();

        gatkCommand.append(java)
                .append(" -Xmx4g -Djava.io.tmpdir=").append(this.tmpDir)
                .append(" -jar ").append(this.gatk_dir).append("/GenomeAnalysisTK.jar")
                .append(" -T SelectVariants")
                .append(" -R ").append(this.genomeFile)
                .append(" -V ").append(inputFile)
                .append(" -o ").append(outputFile)
                .append(" -select 'set == \"Intersection\"'");

        jobGATK.setCommand(gatkCommand.toString());
        jobGATK.setMaxMemory("6000");

        if (!this.queue.isEmpty()) {
            jobGATK.setQueue(this.queue);
        }

        return jobGATK;
    }

    /**
     * This job is for bg-zipping and indexing vcf files with tabix
     *
     * @param inputFile
     * @return
     */
    protected Job prepareVcf(String inputFile) {

        Job jobVcfPrep = wf.createBashJob("prepare_vcf");
        jobVcfPrep.setCommand(getWorkflowBaseDir() + "/dependencies/postprocess_vcf.pl "
                + "--input=" + inputFile + " "
                + "--tabix=" + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/tabix "
                + "--bgzip=" + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/bgzip");
        jobVcfPrep.setMaxMemory("2000");

        if (!this.queue.isEmpty()) {
            jobVcfPrep.setQueue(this.queue);
        }

        return jobVcfPrep;
    }

    /**
     * A utility for making base name of a file
     *
     * @param name
     * @return
     */
    private String makeBasename(final String name) {
        String basename = name.substring(name.lastIndexOf("/") + 1, name.lastIndexOf(".bam"));
        if (basename.contains(".")) {
            basename = basename.substring(0, basename.indexOf("."));
        }
        return basename;
    }

    class InputContainer {

        private String vcfPath;
        private String source;

        public InputContainer() {
            this.vcfPath = null;
            this.source = null;
        }

        public InputContainer(String vcfP, String src) {
            if (!vcfP.contains("vcf.gz")) {
                Log.warn("Files with unconventional extension will be skipped!");
                this.vcfPath = null;
            } else {
                this.vcfPath = vcfP;
                this.source = src;
            }
        }

        public void registerVcf(String vcfP, String src) {
            if (!vcfP.contains("vcf.gz")) {
                Log.warn("Files with unconventional extension will be skipped!");
                return;
            }
            this.vcfPath = vcfP;
            this.source = src;
        }

        /**
         * @return hasVcf
         */
        public boolean hasVcf() {
            return null != this.vcfPath;
        }

        /**
         * @return vcfPath
         */
        public String getVcfPath() {
            return this.vcfPath;
        }

        /**
         * @return source
         */
        public String getSource() {
            return this.source;
        }

    }
}
