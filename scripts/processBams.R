threads <- snakemake@config[['cores']]
library(SpliceWiz)
setSWthreads(threads)
bam_file <- snakemake@input[['aligned_bam']]
ref_path <- snakemake@params[['ref_path']]
processed_bams_path <- snakemake@params[['processed_bams_path']]
aligned_bams_path <- snakemake@params[['aligned_bams_path']]
sample_name = basename(dirname(bam_file))

bams <- findBAMS(aligned_bams_path, level = 1)


processBAM(
    bamfiles = bam_file,
    sample_names = sample_name,
    reference_path = ref_path,
    output_path = processed_bams_path,
    n_threads = threads,
    useOpenMP = TRUE
)
