/**
 * Copyright (C) 2015 Ontario Institute of Cancer Research
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact us:
 *
 * Ontario Institute for Cancer Research MaRS Centre, West Tower 661 University
 * Avenue, Suite 510 Toronto, Ontario, Canada M5G 0A3 Phone: 416-977-7599
 * Toll-free: 1-866-678-6427 www.oicr.on.ca
 *
 */
package ca.on.oicr.pde.deciders;

import java.io.File;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import com.sun.org.apache.xalan.internal.xsltc.compiler.util.Util;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.common.util.maptools.MapTools;
import org.apache.commons.lang3.StringUtils;

public class VariantMergingDecider extends OicrDecider {

    private SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;
    private List<String> variantCallers = null;

    private String output_prefix = "./";
    private String output_dir = "seqware-results";
    private String manual_output = "false";

    private String queue = " ";
    private String sampleName = "NA";
    private String referenceFasta = "";

    private static final String VCF_GZIP_METATYPE  = "application/vcf-4-gzip";
    private static final String VCF_METATYPE       = "text/vcf-4";
    private static final String GVCF_GZIP_METATYPE = "application/g-vcf-gz";
    private static final String CALLER_TOKEN       = "file.variant_caller";
    private static final String TYPE_TOKEN         = "file.variant_type";

    private Map<String, vcfChecker> vcfList;
    private String templateTypeFilter = "";
    private String variantFilter = "";
    private boolean doIntersect;
    private boolean doCollapse;

    public VariantMergingDecider() {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        vcfList        = new HashMap<String, vcfChecker>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("study-name", "Optional*. Set the study name "
                + "will define the scope of the data").withRequiredArg();
        parser.accepts("manual-output", "Optional*. Set the manual output "
                + "either to true or false").withRequiredArg();
        parser.accepts("template-type", "Optional. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type").withRequiredArg();
        parser.accepts("variant-type", "Optional: specify the type of calls (snv, indels) to limit the analysis to this call type").withRequiredArg();
        parser.accepts("output-path", "Optional: the path where the files should be copied to "
                + "after analysis. Corresponds to output-prefix in INI file. Default: ./").withRequiredArg();
        parser.accepts("output-folder", "Optional: the name of the folder to put the output into relative to "
                + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results").withRequiredArg();
        parser.accepts("reference-fasta", "Optional: Reference assembly in fasta format. Default is UCSC hg19_random.fa").withRequiredArg();
        parser.accepts("queue", "Optional: Set the queue (Default: not set)").withRequiredArg();
        parser.accepts("do-intersect", "Optional: Emit only variants called by two SNV callers").withRequiredArg();
        parser.accepts("do-collapse", "Optional: Collapse 4-columns prouced by CombineVariants into two columns (enabled by default)").withRequiredArg();
        parser.accepts("variant-callers", "Optional: Comma-separated list of names for two variant callers that we want to extract data for").withRequiredArg();
        parser.accepts("verbose", "Optional: Enable verbose Logging").withRequiredArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setGroupingStrategy(Header.FILE_SWA);
        ReturnValue val = super.init();

        //Handle .ini file - we accept only memory size allocated to different steps
        if (options.has("ini-file")) {
            File file = new File(options.valueOf("ini-file").toString());
            if (file.exists()) {
                String iniFile = file.getAbsolutePath();
                Map<String, String> iniFileMap = new HashMap<String, String>();
                MapTools.ini2Map(iniFile, iniFileMap);
            } else {
                Log.error("The given INI file does not exist: " + file.getAbsolutePath());
                System.exit(1);
            }

        }

        //Group by sample if no other grouping selected
        if (this.options.has("group-by")) {
            Log.error("group-by parameter passed, but this decider does not allow overriding the default grouping (by FILE_SWA)");
        }

        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        } else {
            this.queue = " ";
        }

        if (this.options.has("template-type")) {
            this.templateTypeFilter = options.valueOf("template-type").toString();
            Log.debug("Setting template type is not necessary, however if set the decider will run the workflow only on this type of data");
        }
        
        if (this.options.has("variant-type")) {
            this.variantFilter = options.valueOf("variant-type").toString();
            Log.debug("Setting variant type is not necessary, however if set the decider will run the workflow only on this type of data");
        }

        if (this.options.has("manual-output")) {
            this.manual_output = options.valueOf("manual_output").toString();
            Log.debug("Setting manual output, default is false and needs to be set only in special cases");
        }

        if (this.options.has("verbose")) {
            Log.setVerbose(true);
        }

        if (this.options.has("output-path")) {
            this.output_prefix = options.valueOf("output-path").toString();
            if (!this.output_prefix.endsWith("/")) {
                this.output_prefix += "/";
            }
        }

        if (this.options.has("output-folder")) {
            this.output_dir = options.valueOf("output-folder").toString();
        }

        if (this.options.has("reference-fasta")) {
            this.referenceFasta = options.valueOf("reference-fasta").toString();
        }

        // Warn about using force-run-all (may not be relevant after 1.0.17 release)
        if (options.has("force-run-all")) {
            Log.stderr("Using --force-run-all WILL BREAK THE LOGIC OF THIS DECIDER, USE AT YOUR OWN RISK");
        }

        // We get NORMAL/TUMOR columns instead of four by default
        if (this.options.has("do-collapse")) {
            this.doCollapse = Boolean.valueOf(options.valueOf("do-collapse").toString());
        } else {
            this.doCollapse = true;
        }

        // We don't want intersect as a default
        if (this.options.has("do-intersect")) {
            this.doIntersect = Boolean.valueOf(options.valueOf("do-intersect").toString());
        } else {
            this.doIntersect = false;
        }

        // Check if we have list of callers, we use it if we won't have caller id metadata 
        if (this.options.has("variant-callers")) {
            this.variantCallers = Arrays.asList(StringUtils.split(options.valueOf("variant-callers").toString(), ","));
            Set<String> inputCallersSet = new HashSet<>(this.variantCallers);

            if (this.variantCallers.size() != inputCallersSet.size()) {
                throw new RuntimeException("Duplicate callers detected in variant-callers");
            }

            if (this.variantCallers.size() != 2) {
                throw new RuntimeException("There must be two callers when merging, variant-callers parameter is invalid");
            }
        }

        // Allows anything defined on the command line to override the defaults here.
        return val;
    }


    @Override
    protected boolean checkFileDetails(FileAttributes fa) {
        Log.debug("CHECK FILE DETAILS:" + fa.getPath());
        String templateType = fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);
        // Filter the data of a different template type if filter is specified
        if (!this.templateTypeFilter.isEmpty() && !this.templateTypeFilter.equalsIgnoreCase(templateType)) {
            Log.stderr("Template type " + templateType + " does not pass filter which is set to " + this.templateTypeFilter);
            return false;
        }

        if (fa.getMetatype().equals(VCF_GZIP_METATYPE) || fa.getMetatype().equals(VCF_METATYPE) || fa.getMetatype().equals(GVCF_GZIP_METATYPE)) {
            BeSmall b = fileSwaToSmall.get(fa.getOtherAttribute(Header.FILE_SWA.getTitle()));

            if (!this.vcfList.containsKey(b.getIusDetails()) || !this.vcfList.get(b.getIusDetails()).hasVcfs() || this.vcfList.get(b.getIusDetails()).getVcfFiles().size() != 2) {
                Log.error("DID NOT FIND suitable vcf files for file [" + b.getPath() + "], SKIPPING");
                return false;
            }
        }
        return true;
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            boolean metatypeOK = false;
            for (int f = 0; f < currentRV.getFiles().size(); f++) {
                try {
                    String type = currentRV.getFiles().get(f).getMetaType();
                    if (type.equals(VCF_GZIP_METATYPE) || type.equals(GVCF_GZIP_METATYPE) || type.equals(VCF_METATYPE)) {

                        metatypeOK = true;
                        FileAttributes fa = new FileAttributes(currentRV, currentRV.getFiles().get(f));

                        // Extract caller info
                        String snvCaller = "";
                        if (null != currentRV.getAttribute(CALLER_TOKEN) && (this.variantCallers == null || this.variantCallers.isEmpty())) {
                            snvCaller = currentRV.getAttribute(CALLER_TOKEN);
                            if (snvCaller.contains(" ")) {
                                continue; // we need to skip these since they are likely overlap files
                            }
                        // As the last resort, try file name for desired caller
                        } else {
                            if (this.variantCallers != null && !this.variantCallers.isEmpty()) {
                                int timeMatched = 0; //Track how many callers the file name matches (to avoid using overlap files)
                                for (String caller : this.variantCallers) {
                                    // Case-agnostic match:
                                    if (Util.baseName(currentRV.getFiles().get(f).getFilePath()).toLowerCase().contains(caller.toLowerCase())) {
                                        snvCaller = caller;
                                        timeMatched++;
                                    }
                                }
                                if (timeMatched != 1) {continue;}
                            }
                        }

                        // We should have extracted caller info by now, if not - skip the file
                        if (null == snvCaller || snvCaller.isEmpty()) {
                            continue;
                        }

                        // Important, the id formed here needs to be exactly as for BeSmall
                        StringBuilder vcfID = new StringBuilder();
                        vcfID.append(fa.getLibrarySample()).append(":")
                                .append(fa.getSequencerRun()).append(":")
                                .append(fa.getLane()).append(":")
                                .append(fa.getBarcode());
                        // Our permissive by-type filter
                        if (null != currentRV.getAttribute(TYPE_TOKEN)) {
                            String varType = currentRV.getAttribute(TYPE_TOKEN);
                            if (!this.variantFilter.isEmpty() && !this.variantFilter.equalsIgnoreCase(varType)) {
                                Log.stderr("Skipping variants of type " + varType);
                                continue;
                            }
                            vcfID.append(":")
                                    .append(currentRV.getAttribute(TYPE_TOKEN));
                        }
                       
                        // Register all vcfs matching this vcfID
                        if (!this.vcfList.containsKey(vcfID.toString())) {
                            this.vcfList.put(vcfID.toString(), new vcfChecker(currentRV.getFiles().get(f).getFilePath(), snvCaller));
                        } else {
                            this.vcfList.get(vcfID.toString()).registerVcf(currentRV.getFiles().get(f).getFilePath(), snvCaller);
                        }
                    }
                } catch (Exception e) {
                    Log.stderr("Error checking a file");
                    continue;
                }
            }

            if (!metatypeOK) {
                continue; // Go to the next value
            }

            BeSmall currentSmall = new BeSmall(currentRV);
            //make sure you only have the most recent single file for each
            //sequencer run + lane + barcode + meta-type
            String fileDeets = currentSmall.getIusDetails();
            Date currentDate = currentSmall.getDate();

            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);           
            
            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                iusDeetsToRV.put(fileDeets, currentRV);
            } //if there is an entry, compare the current value to the 'old' one in
              //the map. if the current date is newer than the 'old' date, replace
              //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);
                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    iusDeetsToRV.put(fileDeets, currentRV);
                }
            }
        }

        //only use those files that entered into the iusDeetsToRV
        //since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        //group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : newValues) {
            String currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).getGroupByAttribute();
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
    }

    @Override
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        BeSmall small = fileSwaToSmall.get(a);
        if (small != null) {
            return small.getGroupByAttribute();
        }
        return attribute;
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {

        StringBuilder groupIds = new StringBuilder();
        String[] filePaths = commaSeparatedFilePaths.split(",");
        StringBuilder tubeId = new StringBuilder();
        StringBuilder groupDescription = new StringBuilder();

        //reset test mode
        if (!this.options.has("test")) {
            this.setTest(false);
        }

        StringBuilder inputVcfs    = new StringBuilder();
        StringBuilder inputSources = new StringBuilder();
        for (BeSmall bs : fileSwaToSmall.values()) {
            if (!bs.getPath().equals(filePaths[0])) {
                continue;
            } else {
                this.sampleName = bs.getSample();
                String[] files   = this.vcfList.get(bs.getIusDetails()).getVcfPaths();
                String[] sources = this.vcfList.get(bs.getIusDetails()).getCallers();
                if (files.length != 2 || sources.length != 2) {
                    Log.stderr("Couldn't get Proper inputs and sources when assembling ini, ABORTING");
                    return null;
                }
                inputVcfs.append(files[0]).append(",").append(files[1]);
                inputSources.append(sources[0]).append(",").append(sources[1]);
                break;
            }
        }

        Map<String, String> iniFileMap = new TreeMap<String, String>();

        iniFileMap.put("input_files",   inputVcfs.toString());
        iniFileMap.put("input_sources", inputSources.toString());
        iniFileMap.put("data_dir", "data");
        iniFileMap.put("output_prefix", this.output_prefix);
        iniFileMap.put("output_dir", this.output_dir);
        iniFileMap.put("identifier", this.sampleName);
        iniFileMap.put("manual_output", this.manual_output);
        iniFileMap.put("collapse_vcf", Boolean.toString(this.doCollapse));
        iniFileMap.put("do_intersect", Boolean.toString(this.doIntersect));
             
        iniFileMap.put("queue", this.queue);
      

        //Note that we can use group_id, group_description and external_name for tumor bams only
        if (null != groupIds && groupIds.length() != 0 && !groupIds.toString().contains("NA")) {
            iniFileMap.put("group_id", groupIds.toString());
        } else {
            iniFileMap.put("group_id", "NA");
        }

        if (null != groupDescription && groupDescription.length() != 0 && !groupIds.toString().contains("NA")) {
            iniFileMap.put("group_id_description", groupDescription.toString());
        } else {
            iniFileMap.put("group_id_description", "NA");
        }

        if (null != tubeId && tubeId.length() != 0 && !groupIds.toString().contains("NA")) {
            iniFileMap.put("external_name", tubeId.toString());
        } else {
            iniFileMap.put("external_name", "NA");
        }

        if (!this.referenceFasta.isEmpty()) {
            iniFileMap.put("reference_fasta", this.referenceFasta);
        }

        return iniFileMap;
    }

    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(VariantMergingDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }

    private class BeSmall {

        private Date   date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String path = null;
        private String sample = null;

        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            StringBuilder iusID = new StringBuilder();

            iusID.append(fa.getLibrarySample()).append(":")
                     .append(fa.getSequencerRun()).append(":")
                     .append(fa.getLane()).append(":")
                     .append(fa.getBarcode());

            if (null != rv.getAttribute(TYPE_TOKEN)) {
                iusID.append(":")
                     .append(rv.getAttribute(TYPE_TOKEN));
            }

            iusDetails = iusID.toString();
            sample = fa.getLibrarySample();

            //if available, use call type information for grouping
            if (null != rv.getAttribute(TYPE_TOKEN)) {
                groupByAttribute = rv.getAttribute(Header.FILE_SWA.getTitle()) + ":" + rv.getAttribute(TYPE_TOKEN);
            } else {
                groupByAttribute = rv.getAttribute(Header.FILE_SWA.getTitle());
            }
            path = rv.getFiles().get(0).getFilePath() + "";
        }

        public Date getDate() {
            return date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getSample() {
            return sample;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }

        public String getIusDetails() {
            return iusDetails;
        }

        public String getPath() {
            return path;
        }

        public void setPath(String path) {
            this.path = path;
        }
    }

    /**
     * This class should have a list of vcf files
     */
    class vcfChecker {

        private Map<String, String> vcfFiles;
        
        /**
         * Public Constructor
         */
        public vcfChecker() {
            this.vcfFiles = new HashMap<String, String>();
        }

        public vcfChecker(String vcfP, String caller) {
            this.vcfFiles = new HashMap<String, String>();
            this.vcfFiles.put(caller, vcfP);
        }

        public void registerVcf(String vcfP, String caller) {
            if (this.vcfFiles.containsKey(caller)) {
                Log.warn("Updating vcf file for " + caller);
            }
            this.vcfFiles.put(caller, vcfP);
        }

        /**
         * @return hasVcf
         */
        public boolean hasVcfs() {
            return !this.vcfFiles.isEmpty();
        }

        /**
         * @return vcfFiles values
         */
        public String[] getVcfPaths() {
            return vcfFiles.values().toArray(new String[0]);
        }
        
        /**
         * @return vcfFiles keySet
         */
        public String[] getCallers() {
            return vcfFiles.keySet().toArray(new String[0]);
        }

        /**
         * @return vcfFiles
         */
        public Map<String, String> getVcfFiles() {
            return this.vcfFiles;
        }
    }

}
