import gzip
# fastq_file = '/home/charlie/projects/splicewiz/input_files/C380_CNS_fastqs/C380_AUB_CNS_rep1_1.fq.gz'
# sam_file = '/home/charlie/projects/filter_bam/filtered_bams/C380_AUB_CNS_rep1_filtered.sam'
# filtered_fastq_file = 'C380_AUB_CNS_rep1_1_filtered.fq'

def get_filtered_reads(sam_file):
	reads = []
	for line in open(sam_file, 'r'):
		line = line.strip().split()
		read = line[0]
		reads.append(read)
	return reads

def write_filtered_fastq(sam_file, fastq_file_1, fastq_file_2, filtered_fastq_1, filtered_fastq_2):
	reads = set(get_filtered_reads(sam_file))
	with gzip.open(fastq_file_1, 'rt') as f, open(filtered_fastq_1, 'w') as out:
		printing = False
		for line in f:
			if line[0] == '@':
				printing = False
				if line.strip().split()[0].strip('@') in reads:
					printing = True
			if printing:
				out.write(line)
	with gzip.open(fastq_file_2, 'rt') as f, open(filtered_fastq_2, 'w') as out:
		printing = False
		for line in f:
			if line[0] == '@':
				printing = False
				if line.strip().split()[0].strip('@') in reads:
					printing = True
			if printing:
				out.write(line)

sam_file = snakemake.input.filtered_sam # type: ignore
fastq_file_1 = snakemake.input.fastq_file_1 # type: ignore
fastq_file_2 = snakemake.input.fastq_file_2 # type: ignore
filtered_fastq_1 = snakemake.output.filtered_fastq_1 # type: ignore
filtered_fastq_2 = snakemake.output.filtered_fastq_2 # type: ignore
write_filtered_fastq(sam_file, fastq_file_1, fastq_file_2, filtered_fastq_1, filtered_fastq_2)