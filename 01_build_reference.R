###############################################################################
#                       build reference variables                             #
###############################################################################
output_folder <- "/home/charlie/splicewiz_output"
cores <- 12
genome_fasta <- "BDGP6.32_fasta/dm_bdgp632.fa"
genome_gtf <- "BDGP6.32_annotations/dm_bdgp632_with_copia.gtf"
ontology_species <- "Drosophila melanogaster"

################################################################################
##                        Build the SpliceWiz reference                       ##
################################################################################
library(SpliceWiz)
setSWthreads(cores)
dir.create(output_folder, showWarnings = FALSE)
ref_path <- file.path(output_folder, gsub(".fa", "", basename(genome_fasta)))
buildRef(
    reference_path = ref_path,
    fasta = genome_fasta, 
    gtf = genome_gtf,
    genome_type = "",
    ontologySpecies = ontology_species
)

################################################################################
##                             Build the STAR reference                       ##
################################################################################
genomeSAindexNbases <- floor(log2(file.info(genome_fasta)$size)/2-1)
STAR_buildRef(
    reference_path = ref_path,
    n_threads = cores,
    additional_args = c("--genomeSAindexNbases", genomeSAindexNbases)
)

