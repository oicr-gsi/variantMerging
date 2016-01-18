/**
 * Copyright (C) 2016 Ontario Institute of Cancer Research
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact us:
 *
 * Ontario Institute for Cancer Research
 * MaRS Centre, West Tower
 * 661 University Avenue, Suite 510
 * Toronto, Ontario, Canada M5G 0A3
 * Phone: 416-977-7599
 * Toll-free: 1-866-678-6427
 * www.oicr.on.ca
 *
 */

package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.SemanticWorkflow;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
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
    private Workflow wf;
    private String tabixVersion;
    private final static String VCF_GZIPPED_METATYPE = "application/vcf-4-gzip";
    private final static String TBI_INDEX_METATYPE   = "application/tbi";
    private boolean manualOutput = false;
    private SqwFile[] inputVcfFiles;
    private ArrayList<String> inputBasenames;
    private ArrayList<String> sources;
    private boolean collapse;

    //References
    private String genomeFile;
    private String genomeDictionary;
    private String queue;
    private String binDir;

    //Ontology-related variables
    private static final String EDAM = "EDAM";
    private static final Map<String, Set<String>> cvTerms;

    //Our CV terms
    static {
        cvTerms = new HashMap<String, Set<String>>();
        cvTerms.put(EDAM, new HashSet<String>(
                          Arrays.asList("VCF", "SNP", "tabix",
                                        "Variant Calling", "Sequence variation analysis")));
    }


    /**
     * Launched from setupDirectory();
     */
    private void VariantMergingWorkflow() {
        identifier = getProperty("identifier");
        java    = getProperty("java");
        dataDir = getOptionalProperty("data_dir", "data");
        tmpDir  = getOptionalProperty("tmp_dir",  "temp");
        binDir  = this.getWorkflowBaseDir() + "/bin/";
        String manual_out = getOptionalProperty("manual_output", "false");
        manualOutput = Boolean.parseBoolean(manual_out);
        wf      = this.getWorkflow();
    }

    @Override
    public Map<String, SqwFile> setupFiles() {

        //References
        this.genomeFile = getProperty("ref_fasta");
        String derivedDictionary = this.genomeFile.substring(0, this.genomeFile.lastIndexOf(".fa")) + ".dict";
        this.genomeDictionary = getOptionalProperty("ref_dictionary", derivedDictionary);

        //GATK
        this.queue    = getOptionalProperty("queue", "");
        this.gatk_dir = getOptionalProperty("gatk_dir", binDir);
        if (!gatk_dir.endsWith("/")) {gatk_dir += "/";}

        //Collapse 4-column final vcf into 2-column vcf
        String collapseColumns = getOptionalProperty("collapse_vcf", "true");
        this.collapse = Boolean.valueOf(collapseColumns);
        
        //Picard
        picard_dir = getOptionalProperty("picard_dir", binDir);
        if (!picard_dir.endsWith("/")) {gatk_dir += "/";}
        
        if (getProperty("tabix_version") == null) {
            Logger.getLogger(VariantMergingWorkflow.class.getName()).log(Level.SEVERE, "tabix_version is not set, we need it to call tabix correctly");
            return (null);
        } else {
            this.tabixVersion = getProperty("tabix_version");
        }

        if (getProperty("java") == null) {
            Logger.getLogger(VariantMergingWorkflow.class.getName()).log(Level.SEVERE, "java is not set, we need it to run GATK since it requires 1.7 java and we cannot rely on the defaul");
            return (null);
        } else {
            this.java = getProperty("java");
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
                this.inputBasenames.add(i, this.makeBasename(this.inputVcfFiles[i].getProvisionedPath()));
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
        StringBuilder infiles = new StringBuilder();
        for (int in = 0; in < this.inputVcfFiles.length; in++) {
            if (in > 0) {infiles.append(",");}
            infiles.append(this.inputVcfFiles[in].getProvisionedPath());
        }
        
        for (int i = 0; i < this.sources.size(); i++) {
            // vet the vcfs
            operationsOnFile = ".vetted.";
            // TODO 
            inputFile = this.inputVcfFiles[i].getProvisionedPath();
            outputFile = this.inputBasenames.get(i) + operationsOnFile;

            Job vetVcf = this.vetVcf(inputFile, this.dataDir + outputFile + "vcf", infiles.toString(), i);

            // update seq dictionary / sort
            operationsOnFile += "sorted.";

            inputFile = outputFile;
            outputFile = this.inputBasenames.get(i) + operationsOnFile;

            Job sortVcf = this.sortVcf(this.dataDir + inputFile + "vcf", this.dataDir + outputFile + "vcf", vetVcf);
            inputsToCombine.put(this.sources.get(i), this.dataDir + outputFile + "vcf");
            upstreamJobs.add(sortVcf);
        }

        /**
         * ==== Combine and Subset (extract intersect) here
         */
        for (int p = 0; p < this.sources.size(); p++) {

            // Combine based on priority
            String operationsOnMergedFile = operationsOnFile + "combined.";
            String priorityOrder = p == 0 ? this.sources.get(0) + "." + this.sources.get(1) + "."
                                          : this.sources.get(1) + "." + this.sources.get(0) + ".";
            outputFile = identifier + operationsOnMergedFile + priorityOrder;
            Job combineVcfs = this.combineVcf(inputsToCombine, this.dataDir + outputFile + "vcf", p);

            for (Job j : upstreamJobs) {
                combineVcfs.addParent(j);
            }
            // TODO : make optional the ability to provision

            // produce intersect (for both overlaps)
            operationsOnMergedFile += "intersect.";
            inputFile = outputFile;
            outputFile = identifier + operationsOnMergedFile + priorityOrder;
            Job subsetVcf = this.subsetVcf(this.dataDir + inputFile + "vcf", this.dataDir + outputFile + "vcf");

            subsetVcf.addParent(combineVcfs);
            
            //post-process vcf (bgzip and index)
            Job postprocessVcf = this.prepareVcf(this.dataDir + outputFile + "vcf");
            postprocessVcf.addParent(subsetVcf);
            
            // Annotate and provision final vcf and its index
            SqwFile finalVcf    = this.createOutputFile(this.dataDir + outputFile + "vcf.gz", VCF_GZIPPED_METATYPE, manualOutput);
            SqwFile finalVcfTbi = this.createOutputFile(this.dataDir + outputFile + "vcf.gz.tbi", TBI_INDEX_METATYPE, manualOutput);
                       
            this.attachCVterms(finalVcf,    EDAM, "VCF,SNP,Variant Calling,Sequence variation analysis");
            this.attachCVterms(finalVcfTbi, EDAM, "tabix,Sequence variation analysis");
            
            //We produce total of four files, vcfs and their indices
            postprocessVcf.addFile(finalVcf);                     
            postprocessVcf.addFile(finalVcfTbi);
            
        }
    }

    //=======================Jobs as functions===================
    /**
     * Vet vcf files - remove tool-specific info from the headers, remove non-canonical contigs
     * TODO need to disambiguate format fields
     * @param inputFile
     * @param outputFile
     * @return
     */
    protected Job vetVcf(String inputFile, String outputFile, String commaSeparatedInputs, int index) {
        
        Job jobVcfPrep = wf.createBashJob("preprocess_vcf");
        jobVcfPrep.setCommand(getWorkflowBaseDir() + "/dependencies/vcf_vetting.pl"
                + " --input "  + inputFile
                + " --infiles "    + commaSeparatedInputs
                + " --index "  + (index + 1) // index comes zero-based, convert to 1-based
                + " --output " + outputFile);

        jobVcfPrep.setMaxMemory("4000");

        if (!this.queue.isEmpty()) {
            jobVcfPrep.setQueue(this.queue);
        }

        return jobVcfPrep;
    }

    /**
     * Sort vcf files using a specified list of contigs
     * @param inputFile
     * @param outputFile
     * @return
     */
    protected Job sortVcf(String inputFile, String outputFile, Job parent) {

        Job jobUpdate = wf.createBashJob("update_dictionary");
        String tempFile = inputFile.replace(".vcf", ".updated.vcf");

        jobUpdate.setCommand(this.java + " -Xmx" + this.getProperty("picard_memory") + "M"
                + " -jar " + this.picard_dir + "picard.jar UpdateVcfSequenceDictionary "
                + " I=" + inputFile
                + " O=" + tempFile
                + " SEQUENCE_DICTIONARY=" + this.genomeDictionary
                + " CREATE_INDEX=false");

        Integer mbRam = Integer.valueOf(this.getProperty("picard_memory"));
        jobUpdate.setMaxMemory(Integer.toString(mbRam + 2000));
        jobUpdate.addParent(parent);

        Job jobSort = wf.createBashJob("sort_variants");

        jobSort.setCommand(this.java + " -Xmx" + this.getProperty("picard_memory") + "M"
                + " -jar " + this.picard_dir + "picard.jar SortVcf"
                + " R=" + this.genomeFile
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
     * Combine the vcf files into one according to priorities
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

        gatkCommand.append(this.java)
                .append(" -Xmx4g -Djava.io.tmpdir=").append(this.tmpDir)
                .append(" -jar ").append(this.gatk_dir).append("GenomeAnalysisTK.jar")
                .append(" -T CombineVariants ")
                .append(" --variant:").append(this.sources.get(priorityIndex)).append(" ").append(inputs.get(this.sources.get(priorityIndex)))
                .append(" --variant:").append(this.sources.get(otherIndex)).append(" ").append(inputs.get(this.sources.get(otherIndex)))
                .append(" -R ").append(this.genomeFile)
                .append(" -o ").append(outputFile)
                .append(" -genotypeMergeOptions UNIQUIFY"); // RIORITIZE")
                //.append(" -priority ").append(this.sources.get(priorityIndex)).append(",").append(this.sources.get(otherIndex));

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

        gatkCommand.append(this.java)
                .append(" -Xmx4g -Djava.io.tmpdir=").append(this.tmpDir)
                .append(" -jar ").append(this.gatk_dir).append("GenomeAnalysisTK.jar")
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
        //jobVcfPrep.setCommand
        StringBuilder postprocessCommand = new StringBuilder();
        postprocessCommand.append(getWorkflowBaseDir())
                .append("/dependencies/postprocess_vcf.pl")
                .append(" --input ").append(inputFile)
                .append(" --tabix ").append(getWorkflowBaseDir()).append("/bin/tabix-").append(this.tabixVersion).append("/tabix")
                .append(" --bgzip ").append(getWorkflowBaseDir()).append("/bin/tabix-").append(this.tabixVersion).append("/bgzip");
        if (this.collapse) {
            postprocessCommand.append(" --collapse ");
        }
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
        return name.substring(name.lastIndexOf("/") + 1, name.lastIndexOf(".vcf"));
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

}
