library(SpliceWiz)
ref_path <- snakemake@params[['ref_path']]
cores <- snakemake@config[['cores']]
setSWthreads(cores)

################################################################################
##                        Build the SpliceWiz reference                       ##
################################################################################

buildRef(
    reference_path = ref_path,
    fasta = "", 
    gtf = "",
    genome_type = "",
)

