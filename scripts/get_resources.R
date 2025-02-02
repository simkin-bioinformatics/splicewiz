cores <- snakemake@config[['cores']]
library(SpliceWiz)
setSWthreads(cores)
ref_path <- snakemake@params[['ref_path']]
genome_fasta <- snakemake@input[['genome_fasta']]
genome_gtf <- snakemake@input[['genome_gtf']]

getResources(
    reference_path = ref_path,
    fasta = genome_fasta,
    gtf = genome_gtf
)
