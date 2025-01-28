import os
aligned_bams_folder = "/home/charlie/big_data/splicewiz_output/dm_bdgp632/C380_CNS_alignment/aligned_bams"
fastq_folder = "/home/charlie/Dropbox/splicewiz/splicewiz_fastqs/C380_CNS"

bed_file = "/home/charlie/big_data/genomes_and_TE/BDGP6.32_annotations/dm_bdgp632_copia_regions.bed"


output_folder = f"{fastq_folder}_filtered"
samples = set({})
for file in os.listdir(fastq_folder):
	file = file.replace("_1.fq.gz","").replace("_2.fq.gz", "")
	# file = file.strip("_2.fq.gz")
	samples.add(file)
samples = list(samples)
samples.sort()

rule all:
	input:
		expand(os.path.join(output_folder, "{sample}_1.fq.gz"), sample = samples),
		expand(os.path.join(output_folder, "{sample}_2.fq.gz"), sample = samples),
rule sort_bams:
	input: 
		os.path.join(aligned_bams_folder, "{sample}", "Aligned.out.bam")
	output: 
		os.path.join(output_folder, "sorted_bams", "{sample}_sorted.bam")
	shell: 
		"samtools sort -o {output} {input}"

rule index_bams:
	input:
		os.path.join(output_folder, "sorted_bams", "{sample}_sorted.bam")
	output:
		os.path.join(output_folder, "sorted_bams", "{sample}_sorted.bam.bai")
	shell:
		"samtools index {input}"

rule filter_bams:
	input:
		indexes = os.path.join(output_folder, "sorted_bams", "{sample}_sorted.bam.bai"),
		bams = os.path.join(output_folder, "sorted_bams", "{sample}_sorted.bam"),
		bed_file = bed_file
	output:
		filtered_sam = os.path.join(output_folder, "filtered_bams", "{sample}_filtered.sam")
	shell:
		"samtools view -o {output} {input.bams} -L {input.bed_file}"

rule filter_fastqs:
	input:
		fastq_file_1 = os.path.join(fastq_folder, "{sample}_1.fq.gz"),
		fastq_file_2 = os.path.join(fastq_folder, "{sample}_2.fq.gz"),
		filtered_sam = os.path.join(output_folder, "filtered_bams", "{sample}_filtered.sam")
	output:
		filtered_fastq_1 = os.path.join(output_folder, "{sample}_1.fq"),
		filtered_fastq_2 = os.path.join(output_folder, "{sample}_2.fq"),
	script:
		"filter_fastqs.py"

rule compress_filtered_fastqs:
	input:
		filtered_fastq_1 = os.path.join(output_folder, "{sample}_1.fq"),
		filtered_fastq_2 = os.path.join(output_folder, "{sample}_2.fq")
	output:
		compressed_fastq_1 = os.path.join(output_folder, "{sample}_1.fq.gz"),
		compressed_fastq_2 = os.path.join(output_folder, "{sample}_2.fq.gz"),
	shell:
		"pigz {input}"
