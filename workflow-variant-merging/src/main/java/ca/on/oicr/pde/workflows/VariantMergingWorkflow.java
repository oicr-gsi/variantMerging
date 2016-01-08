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

    // TODO review these parameters
    private final List<String> inputVcfFiles = new LinkedList<>();
    private String finalOutput;

    private String dataDir, tmpDir;
    private String java;
    private boolean manualOutput = false;
    private String picard_dir;
    private String gatk_dir;
    private String identifier;
    private String sourceA;
    private String sourceB;
    private String alignerName;
    private Workflow wf;
    private PicardTools picard;
    private static final String REMOVE_DUPLICATES = "REMOVE_DUPLICATES=";
    private final static String VCF_GZIPPED_METATYPE = "application/vcf-4-gzip";
    private final static String VCF_METATYPE = "application/vcf-4";

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
        cvTerms.put(EDAM, new HashSet<String>(Arrays.asList("VCF", "plain text format (unformatted)",
                "Read pre-processing", "Sequence alignment refinement",
                "Sequence alignment metadata")));
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
        picard = new PicardTools(wf);
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

        if (getProperty("gatk_java") == null) {
            Logger.getLogger(VariantMergingWorkflow.class.getName()).log(Level.SEVERE, "gatk_java is not set, we need it to run GATK since it requires 1.7 java and we cannot rely on the defaul");
            return (null);
        } else {
            this.java = getProperty("gatk_java");
        }
        List<String> inputFilesList = Arrays.asList(StringUtils.split(getProperty("input_files"), ","));
        Set<String> inputFilesSet = new HashSet<>(inputFilesList);

        if (inputFilesList.size() != inputFilesSet.size()) {
            throw new RuntimeException("Duplicate files detected in input_files");
        }

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
        String operationsOnFile = identifier + ".";
        String outputFile = "";
        String inputFile = outputFile;
        List<Job> upstreamJobs = new ArrayList<Job>();

        //=======TODO
        // vet the vcfs
        operationsOnFile += "vetted.";

        // update seq dictionary / sort
        operationsOnFile += "sorted.";

        // combine
        operationsOnFile += "combined.";

        // produce intersect
        operationsOnFile += "selected.";

        // Annotate and provision final bam and its index
        SqwFile finalVcf = this.createOutputFile(this.finalOutput + "vcf.gz", VCF_GZIPPED_METATYPE, manualOutput);
        if (!this.alignerName.isEmpty()) {
            finalVcf.getAnnotations().put("aligner", this.alignerName);
        }

        this.attachCVterms(finalVcf, EDAM, "BAM,Sequence alignment refinement");
       // jobMergeFinal.addFile(finalVcf);

        //Job jobIdx2 = this.indexBamJob(this.finalOutput);
        //jobIdx2.addParent(jobMergeFinal);
    }

    //=======================Jobs as functions===================
    /**
     *
     * @param inputFile
     * @param outputFile
     * @return
     */
    protected Job vetVcf(String inputFile, String outputFile) {
        Job jobVcfPrep   = wf.createBashJob("preprocess_vcf");
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
        
        Job jobUpdate   = wf.createBashJob("update_dictionary");
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
        return jobSort;
    }

    /**
     * combine the vcf files into one according to priorities
     *
     * @param inputs
     * @param outputFile
     * @return
     */
    protected Job combineVcf(Map<String, String> inputs, String outputFile) {
        Job jobGATK = this.wf.createBashJob("combine_calls");
        StringBuilder gatkCommand = new StringBuilder();

        gatkCommand.append(java)
                .append(" -Xmx4g -Djava.io.tmpdir=").append(this.tmpDir)
                .append(" -jar ").append(this.gatk_dir).append("/GenomeAnalysisTK.jar")
                .append(" -T CombineVariants ")
                .append(" --variant:").append(this.sourceA).append(" ").append(inputs.get(this.sourceA))
                .append(" --variant:").append(this.sourceB).append(" ").append(inputs.get(this.sourceB))
                .append(" -R ").append(this.genomeFile)
                .append(" -o ").append(outputFile)
                .append(" -genotypeMergeOptions PRIORITIZE")
                .append(" -priority ").append(this.sourceA).append(",").append(this.sourceB);

        jobGATK.setCommand(gatkCommand.toString());
        jobGATK.setMaxMemory("6000");
        return jobGATK;
    }

    /**
     * Produce intersect subset of the calls
     *
     * @param inputs
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
        return jobGATK;
    }
}
