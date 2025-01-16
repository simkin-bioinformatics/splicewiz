################################################################################
##                          Experiment Visualization                          ##
################################################################################

# load splicewiz library and config file
library(SpliceWiz)
config <- yaml::yaml.load_file("config.yml")
ref_path = file.path(config$output_folder, gsub(".fa", "", basename(config$genome_fasta)))


# set the input and output folders
analysis_path <- file.path(
    ref_path, 
    paste("analysis", basename(config$fastq_folder_path), sep='_')
)
vis_path <- file.path(analysis_path, "visualization")

# load the objects from analysis
se <- readRDS(file = file.path(vis_path, "se_high_confidence.RDS"))
res_edgeR <- readRDS(file = file.path(vis_path, "res_edgeR_high_confidence.RDS"))

# create coverage plot
dataObjCopia <- getCoverageData(se, Gene = config$gene, tracks = colnames(se))
plotObjCopia <- getPlotObject(dataObjCopia, Event = res_edgeR$EventName[config$row])
coverage_plot_name <- paste(config$gene, config$row, sep = "_")
pdf(
    file.path(vis_path,paste(coverage_plot_name, ".pdf", sep="")),
    width=10,
    height=20
)
plotView(
    plotObjCopia,
    centerByEvent = TRUE, # whether the plot should be centered at the `Event`
    trackList = list(1,2,3,4,5,6),
    plotJunctions = TRUE
)
dev.off()
