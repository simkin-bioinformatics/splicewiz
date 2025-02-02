library(SpliceWiz)
library(glue)
# chroms: full, Y, X, 3R, 3L, 2R, 2L
###############################################################################
#                       build reference variables                             #
###############################################################################
sample <- "C380_CNS"
fastq_type <- "full"
for (chrom in c('full', 'Y', 'X', '3R', '3L', '2R', '2L')){
# fastq_type <- "full"
# for (chrom in c('full')){
genome_fasta <- glue("/home/charlie/projects/reference/genome_fastas/BDGP6.46.{chrom}.fa")
genome_gtf <- "/home/charlie/projects/reference/genome_annotations/BDGP6.46.full_with_copia.gtf"

fastq_folder <- glue("/home/charlie/projects/reference/splicewiz_fastqs/{sample}_{fastq_type}")
fastq_suffix <- ".fq.gz"
is_paired <- TRUE
################################################################################
genome_name = gsub(".fa", "",basename(genome_fasta))
sample_name = basename(fastq_folder)
ref_path <- file.path("output", glue('reference_{genome_name}'))
alignment_path <- file.path()
aligned_bams_path <- file.path(ref_path, glue('aligned_{sample_name}'), 'aligned_bams')
processed_bams_path <- gsub('aligned_bams', 'processed_bams', aligned_bams_path)
nxtse_path <- gsub('aligned_bams', 'nxtse', aligned_bams_path)
dir.create(ref_path, recursive = TRUE, showWarnings=FALSE)
dir.create(aligned_bams_path, recursive = TRUE, showWarnings=FALSE)
dir.create(processed_bams_path, recursive = TRUE, showWarnings=FALSE)
dir.create(nxtse_path, recursive = TRUE, showWarnings=FALSE)

################################################################################
##                        Build the SpliceWiz reference                       ##
################################################################################

try(buildRef(
    reference_path = ref_path,
    fasta = genome_fasta, 
    gtf = genome_gtf,
    genome_type = "",
))

################################################################################
##                             Build the STAR reference                       ##
################################################################################
genomeSAindexNbases <- floor(log2(file.info(genome_fasta)$size)/2-1)
try(STAR_buildRef(
    reference_path = ref_path,
    n_threads = 8,
    additional_args = c("--genomeSAindexNbases", genomeSAindexNbases)
))

################################################################################
##                               STAR alignment                               ##
################################################################################

fastq_files <- findFASTQ(
    sample_path = fastq_folder, 
    paired = is_paired,
    fastq_suffix = fastq_suffix,
    level = 0
)

STAR_alignExperiment(
    Experiment = fastq_files,
    STAR_ref_path = file.path(ref_path, 'STAR'),
    BAM_output_path = aligned_bams_path,
    n_threads = 8,
    trim_adaptor = "",
    additional_args = c("--outFilterMultimapNmax", "40"), #in case a read with copia maps to many places
    two_pass = FALSE
)

################################################################################
##                                   Process Bams                             ##
################################################################################

bams <- findBAMS(aligned_bams_path, level = 1)

# process aligned bams
processBAM(
    bamfiles = bams$path,
    sample_names = bams$sample,
    reference_path = ref_path,
    output_path = processed_bams_path,
    n_threads = 8,
    useOpenMP = TRUE
)

################################################################################
##                           Collate the experiment                           ##
################################################################################

collateData(
    Experiment = findSpliceWizOutput(processed_bams_path),
    reference_path = ref_path,
    output_path = nxtse_path,
    novelSplicing = FALSE
)


}
