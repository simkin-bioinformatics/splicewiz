################################################################################
##                        Build the splicewiz reference                       ##
################################################################################

library(SpliceWiz)
config <- yaml::yaml.load_file("config.yml")
ref_path = file.path(config$output_folder, gsub(".fa", "", basename(config$genome_fasta)))
if (!dir.exists(config$output_folder)) {dir.create(config$output_folder)}
setSWthreads(config$cores)
buildRef(
    reference_path = ref_path,
    fasta = config$genome_fasta, 
    gtf = config$genome_gtf,
    genome_type = "",
    ontologySpecies = config$ontology_species,
)

################################################################################
##                             Build the STAR reference                       ##
################################################################################
library(SpliceWiz)
config <- yaml::yaml.load_file("config.yml")
ref_path = file.path(config$output_folder, gsub(".fa", "", basename(config$genome_fasta)))
setSWthreads(config$cores)
STAR_buildRef(
    reference_path = ref_path,
    n_threads = config$cores,
    additional_args = c("--genomeSAindexNbases", config$genomeSAindexNbases)
)

