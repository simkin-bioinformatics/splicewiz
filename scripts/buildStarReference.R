library(SpliceWiz)
ref_path <- snakemake@params[['ref_path']]
cores <- snakemake@config[['cores']]
setSWthreads(cores)
genome_fasta <- snakemake@config[['genome_fasta']]

genomeSAindexNbases <- floor(log2(file.info(genome_fasta)$size)/2-1)
STAR_buildRef(
    reference_path = ref_path,
    n_threads = 8,
    additional_args = c("--genomeSAindexNbases", genomeSAindexNbases),
    overwrite = TRUE
)
