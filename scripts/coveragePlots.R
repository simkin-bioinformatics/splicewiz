library(SpliceWiz)
library(glue)
setSWthreads(snakemake@config['cores'])
###############################################################################
out_name <- snakemake@params[['out_name']]
nxtse_path <- snakemake@params[['nxtse_path']]
output_folder <- snakemake@config[['output_folder']]
deltaPSI_cutoff <- snakemake@config[['deltaPSI_cutoff']]
search_term <- snakemake@config[['search_term']]
alignment_path <- snakemake@params[['alignment_path']]
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
write.csv(res_edgeR_searched, file.path(alignment_path, glue('{out_name}_resEdgR.csv')))
for (event in res_edgeR_searched$EventName) {
    dataObj <- getCoverageData(se, Event=event, tracks = colnames(se))
    plotObj <- getPlotObject(dataObj, Event = event)
    gene = strsplit(event, split = "/")[[1]][1]
    coverage_plot_name <- glue("{out_name}_coverage_{gene}.pdf")
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
