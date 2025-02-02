library(SpliceWiz)
setSWthreads(snakemake@config[['cores']])
fastq_folder <- snakemake@config[['fastq_folder']]
is_paired <- snakemake@config[['is_paired']]
fastq_suffix <- snakemake@config[['fastq_suffix']]
ref_path <- snakemake@params[['ref_path']]
aligned_bams_path <- snakemake@params[['aligned_bams_path']]

fastq_1 <- snakemake@input[['fastq_1']]
fastq_2 <- snakemake@input[['fastq_2']]


STAR_alignReads(
    fastq_1 = fastq_1, fastq_2 = fastq_2,
    STAR_ref_path = file.path(ref_path, "STAR"),
    BAM_output_path = aligned_bams_path,
    n_threads = snakemake@config[['cores']],
    trim_adaptor = "",
    additional_args = c("--outFilterMultimapNmax", "40"), #in case a read with copia maps to many places
    two_pass = FALSE
)
