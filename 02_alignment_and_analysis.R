###############################################################################
#    choose an output folder (reference folder will be created inside)        #
###############################################################################
ref_path <- "output/BDGP6.32.full"
cores <- 30
fastq_folder <- "/home/charlie/projects/reference/splicewiz_fastqs/C380_CNS"
sample_name <- "C380_CNS"
fastq_suffix <- ".fq.gz" # e.g. .fq or .fq.gz
is_paired <- TRUE # TRUE or FALSE


deltaPSI_cutoff <- 0.00
search_term <- "copia"
################################################################################
##                               STAR alignment                               ##
################################################################################
library(SpliceWiz)
library(glue)
setSWthreads(cores)
sample <- basename(fastq_folder)
alignment_path <- file.path(ref_path, glue("{sample}_alignment"))
aligned_bams_path <- file.path(alignment_path, "aligned_bams")
dir.create(aligned_bams_path, recursive = TRUE, showWarnings = FALSE)

fastq_files <- findFASTQ(
    sample_path = fastq_folder, 
    paired = is_paired,
    fastq_suffix = fastq_suffix,
    level = 0
)

STAR_alignExperiment(
    Experiment = fastq_files,
    STAR_ref_path = file.path(ref_path, "STAR"),
    BAM_output_path = aligned_bams_path,
    n_threads = cores,
    trim_adaptor = "",
    additional_args = c("--outFilterMultimapNmax", "40"), #in case a read with copia maps to many places
    two_pass = FALSE
)

################################################################################
##                               Process BAM files                            ##
################################################################################
processed_bams_path <- file.path(alignment_path, "processed_bams")
dir.create(processed_bams_path, recursive = TRUE, showWarnings = FALSE)

bams <- findBAMS(aligned_bams_path, level = 1)

# process aligned bams
processBAM(
    bamfiles = bams$path,
    sample_names = bams$sample,
    reference_path = ref_path,
    output_path = processed_bams_path,
    n_threads = cores,
    useOpenMP = TRUE
)

################################################################################
##                           Collate the experiment                           ##
################################################################################

nxtse_path <- file.path(alignment_path, "nxtse")
dir.create(nxtse_path, recursive = TRUE, showWarnings = FALSE)

# create the nxtse files
collateData(
    Experiment = findSpliceWizOutput(processed_bams_path),
    reference_path = ref_path,
    output_path = nxtse_path,
    novelSplicing = FALSE
)

################################################################################
##                          Experiment analysis                               ##
################################################################################
# import the experiment
se <- makeSE(nxtse_path)

# annotate the experiment
colData(se)$condition <- rep(c("AUB", "CS"), each = 3)

# Apply default filters
se_filtered <- se[applyFilters(se),]

# Perform Differential Analysis
res_edgeR <- ASE_edgeR(
    se = se_filtered,
    test_factor = "condition",
    test_nom = "AUB",
    test_denom = "CS"
)

################################################################################
##                            View Coverage Plot                              ##
################################################################################
res_edgeR_filtered <- res_edgeR[res_edgeR$abs_deltaPSI > deltaPSI_cutoff,]
res_edgeR_searched <- res_edgeR_filtered[grepl(search_term, res_edgeR_filtered$EventName, ignore.case = TRUE)]
write.csv(res_edgeR_filtered, file.path(alignment_path, glue('{sample_name}_res_edgeR.csv')))
write.csv(res_edgeR_searched, file.path(alignment_path, glue('{sample_name}_res_edgeR_searched.csv')))
for (event in res_edgeR_searched$EventName) {
    dataObj <- getCoverageData(se, Event=event, tracks = colnames(se))
    plotObj <- getPlotObject(dataObj, Event = event)
    gene = strsplit(event, split = "/")[[1]][1]
    coverage_plot_name <- glue("{basename(ref_path)}_{sample_name}_{gene}.pdf")
        pdf(
            file.path(alignment_path,coverage_plot_name),
            width=10,
            height=20
        )
        try(print(plotView(
            plotObj,
            centerByEvent = TRUE, # whether the plot should be centered at the `Event`
            trackList = list(1,2,3,4,5,6),
            plotJunctions = TRUE
        )))
        dev.off()
}

