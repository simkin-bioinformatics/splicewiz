library(SpliceWiz)
processed_bams_path <- snakemake@params[['processed_bams_path']]
ref_path <- snakemake@params[['ref_path']]
nxtse_path <- snakemake@params[['nxtse_path']]
collateData(
    Experiment = findSpliceWizOutput(processed_bams_path),
    reference_path = ref_path,
    output_path = nxtse_path,
    novelSplicing = FALSE
)
