###############################################################################
#                       build reference variables                             #
###############################################################################
output_folder <- "output"
cores <- 30
genome_fasta <- "/home/charlie/projects/reference/genome_fastas/BDGP6.32.full.fa"
genome_gtf <- "/home/charlie/projects/reference/genome_annotations/BDGP6.32.full_plus_copia.gtf"
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

