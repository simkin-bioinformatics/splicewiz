################################################################################
##                               STAR alignment                               ##
################################################################################
# load splicewiz library and config file
library(SpliceWiz)
config <- yaml::yaml.load_file("config.yml")
ref_path = file.path(config$output_folder, gsub(".fa", "", basename(config$genome_fasta)))
setSWthreads(config$cores)

# set paths to output folders and create them if they don't exit
analysis_path <- file.path(
    ref_path, 
    paste("analysis", basename(config$fastq_folder_path), sep='_')
)
aligned_bams_path = file.path(analysis_path, "aligned_bams")
if (!dir.exists(analysis_path)) {dir.create(analysis_path)}
if (!dir.exists(aligned_bams_path)) {dir.create(aligned_bams_path)}

# load fastq files
fastq_files <- findFASTQ(
    sample_path = config$fastq_folder_path, 
    paired = config$is_paired,
    fastq_suffix = config$fastq_suffix,
    level = 1
)

# perform alignent
STAR_alignExperiment(
    Experiment = fastq_files,
    STAR_ref_path = file.path(ref_path, "STAR"),
    BAM_output_path = aligned_bams_path,
    n_threads = config$cores,
    trim_adaptor = "",
    additional_args = c("--outFilterMultimapNmax", "40"), #in case a read with copia maps to many places
    two_pass = FALSE
)

################################################################################
##                               Process BAM files                            ##
################################################################################

# load splicewiz library and config file
library(SpliceWiz)
config <- yaml::yaml.load_file("config.yml")
ref_path = file.path(config$output_folder, gsub(".fa", "", basename(config$genome_fasta)))
setSWthreads(config$cores)

# set the output folders and create them if needed
analysis_path <- file.path(
    ref_path, 
    paste("analysis", basename(config$fastq_folder_path), sep='_')
)
aligned_bams_path <- file.path(analysis_path, "aligned_bams")
processed_bams_path <- file.path(analysis_path, "processed_bams")
bams <- findBAMS(aligned_bams_path, level = 1)
if (!dir.exists(processed_bams_path)) {dir.create(processed_bams_path)}

# process aligned bams
processBAM(
    bamfiles = bams$path,
    sample_names = bams$sample,
    reference_path = ref_path,
    output_path = processed_bams_path,
    n_threads = config$cores,
    useOpenMP = TRUE
)

################################################################################
##                           Collate the experiment                           ##
################################################################################

# load splicewiz library and config file
library(SpliceWiz)
config <- yaml::yaml.load_file("config.yml")
ref_path = file.path(config$output_folder, gsub(".fa", "", basename(config$genome_fasta)))
setSWthreads(config$cores)

# set the output folders and create them if needed
analysis_path <- file.path(
    ref_path, 
    paste("analysis", basename(config$fastq_folder_path), sep='_')
)
processed_bams_path <- file.path(analysis_path, "processed_bams")
nxtse_path <- file.path(analysis_path, "nxtse")
if (!dir.exists(nxtse_path)) {dir.create(nxtse_path)}

# create the nxtse files
collateData(
    Experiment = findSpliceWizOutput(processed_bams_path),
    reference_path = ref_path,
    output_path = nxtse_path,
    novelSplicing = FALSE
)
