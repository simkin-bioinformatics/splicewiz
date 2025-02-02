configfile: "config.yaml"

import os
genome_name = os.path.basename(config['genome_fasta'].replace('.fa',''))
ref_path = os.path.join(config['output_folder'], f"{os.path.basename(config['genome_fasta'].replace('.fa', ''))}_reference")
sample_name = os.path.basename(config['fastq_folder'])
alignment_path= os.path.join(config['output_folder'], f"{sample_name}_{genome_name}_alignment")
aligned_bams_path = os.path.join(alignment_path, 'aligned_bams')
processed_bams_path = os.path.join(alignment_path, 'processed_bams')
nxtse_path = os.path.join(alignment_path, 'nxtse')

# os.makedirs(os.path.join(ref_path, 'STAR'), exist_ok=True)

rule all:
  input:
    plots_done = os.path.join(ref_path, 'plots_done.txt')
   
rule getResources:
  input:
    genome_fasta = config['genome_fasta'],
    genome_gtf = config['genome_gtf']
  params:
    ref_path = ref_path
  output:
    copied_genome = os.path.join(ref_path, 'resource', 'genome.2bit'),
    copied_gtf = os.path.join(ref_path, 'resource', 'transcripts.gtf.gz')
  script:
    "scripts/get_resources.R"
    
rule SpliceRef:
  input:
    copied_genome = os.path.join(ref_path, 'resource', 'genome.2bit'),
    copied_gtf = os.path.join(ref_path, 'resource', 'transcripts.gtf.gz')
  params:
    ref_path = ref_path
  output:
    splicewiz_ref = os.path.join(ref_path, 'SpliceWiz.ref.gz')
  script:
    "scripts/buildSpliceReference.R"

rule StarRef:
  input:
    copied_genome = os.path.join(ref_path, 'resource', 'genome.2bit'),
    copied_gtf = os.path.join(ref_path, 'resource', 'transcripts.gtf.gz')
  params:
    ref_path = ref_path
  output:
    star_done = touch(os.path.join(ref_path, "STAR", 'done.txt'))
  script:
    "scripts/buildStarReference.R"

rule StarAlignment:
  input:
    star_done = os.path.join(ref_path, "STAR", 'done.txt'),
    fastq_1 = os.path.join(config['fastq_folder'], "{sample}_1.fq.gz"),
    fastq_2 = os.path.join(config['fastq_folder'], "{sample}_2.fq.gz"),
  params:
    aligned_bams_path = os.path.join(aligned_bams_path, '{sample}'),
    ref_path = ref_path
  output:
    aligned_bam = os.path.join(aligned_bams_path, "{sample}", "Aligned.out.bam")
  script:
    "scripts/starAlignment.R"

rule ProcessBams:
  input:
    aligned_bam = os.path.join(aligned_bams_path, "{sample}", "Aligned.out.bam"),
    splicewiz_ref = os.path.join(ref_path, 'SpliceWiz.ref.gz')
  output:
    processed = os.path.join(processed_bams_path, "{sample}.cov")
  params:
    processed_bams_path = processed_bams_path,
    aligned_bams_path = aligned_bams_path,
    ref_path = ref_path,
    sample_name = '{sample}'
  script:
    "scripts/processBams.R"

rule CollateData:
  input:
    processed = expand(os.path.join(processed_bams_path, "{sample}.cov"), sample = config['samples'])
  output:
    filtered_introns = os.path.join(nxtse_path, 'filteredIntrons.Rds')
  params:
    ref_path = ref_path,
    nxtse_path = nxtse_path,
    processed_bams_path = processed_bams_path
  script:
    "scripts/collateData.R"

rule CoveragePlots:
  input:
    filtered_introns = os.path.join(nxtse_path, 'filteredIntrons.Rds')
  output:
    plots_done = temp(touch(os.path.join(ref_path, 'plots_done.txt')))
  params:
    nxtse_path = nxtse_path,
    out_name = f"{sample_name}_{genome_name}",
    ref_path = ref_path,
    alignment_path = alignment_path
  script:
    "scripts/coveragePlots.R"
