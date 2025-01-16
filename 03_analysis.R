################################################################################
##                          Experiment analysis                               ##
################################################################################

# load splicewiz library and config file
library(SpliceWiz)
config <- yaml::yaml.load_file("config.yml")
ref_path = file.path(config$output_folder, gsub(".fa", "", basename(config$genome_fasta)))
setSWthreads(config$cores)

# set the input and output folders
analysis_path <- file.path(
    ref_path, 
    paste("analysis", basename(config$fastq_folder_path), sep='_')
)
nxtse_path <- file.path(analysis_path, "nxtse")
vis_path <- file.path(analysis_path, "visualization")
if (!dir.exists(vis_path)) {dir.create(vis_path)}

# import the experiment
se_unfiltered <- makeSE(nxtse_path)

# annotate the experiment
colData(se_unfiltered)$condition <- rep(c("AUB", "CS"), each = 3)

# Filter high confidence events                    
se_high_confidence <- se_unfiltered[applyFilters(se_unfiltered),]

# Perform Differential Analysis
res_edgeR_high_confidence <- ASE_edgeR(
    se = se_high_confidence,
    test_factor = "condition",
    test_nom = "AUB",
    test_denom = "CS"
)

res_edgeR_unfiltered <- ASE_edgeR(
    se = se_unfiltered,
    test_factor = "condition",
    test_nom = "AUB",
    test_denom = "CS"
)

# write csv of res_edgeR and save objects needed for visualization
write.csv(res_edgeR_unfiltered, file.path(vis_path, 'res_edgeR_unfiltered.csv'))
write.csv(res_edgeR_high_confidence, file.path(vis_path, 'res_edgeR_high_confidence.csv'))
saveRDS(se_high_confidence, file = file.path(vis_path, "se_high_confidence.RDS"))
saveRDS(res_edgeR_high_confidence, file = file.path(vis_path, "res_edgeR_high_confidence.RDS"))
