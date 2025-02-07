// Load configuration from config.yaml
params.config = "config.yaml"
config = file(params.config)

// Define paths and variables
genome_name = file(config.genome_fasta).getName().replace(".fa", "").replace(".gz", "")
sample_name = file(config.fastq_folder).getName()
ref_path = file("${config.output_folder}/${genome_name}_-_reference")
alignment_path = file("${config.output_folder}/${sample_name}_-_{genome_name}_-_alignment")
aligned_bams_path = file("${alignment_path}/aligned_bams")
processed_bams_path = file("${alignment_path}/processed_bams")
nxtse_path = file("${alignment_path}/nxtse")

// Define processes

process getResources {
    input:
    path genome_fasta
    path genome_gtf

    output:
    path "${ref_path}/resource/genome.2bit", emit: copied_genome
    path "${ref_path}/resource/transcripts.gtf.gz", emit: copied_gtf

    script:
    """
    mkdir -p ${ref_path}/resource
    scripts/get_resources.R
    """
}

process SpliceRef {
    input:
    path copied_genome
    path copied_gtf

    output:
    path "${ref_path}/SpliceWiz.ref.gz", emit: splicewiz_ref

    script:
    """
    scripts/buildSpliceReference.R
    """
}

process StarRef {
    input:
    path copied_genome
    path copied_gtf

    output:
    path "${ref_path}/STAR/done.txt", emit: star_done

    script:
    """
    mkdir -p ${ref_path}/STAR
    scripts/buildStarReference.R
    touch ${ref_path}/STAR/done.txt
    """
}

process StarAlignment {
    input:
    path star_done
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    path "${aligned_bams_path}/${sample}/Aligned.out.bam", emit: aligned_bam

    script:
    """
    mkdir -p ${aligned_bams_path}/${sample}
    scripts/starAlignment.R
    """
}

process ProcessBams {
    input:
    path aligned_bam
    path splicewiz_ref

    output:
    path "${processed_bams_path}/${sample}.cov", emit: processed

    script:
    """
    scripts/processBams.R
    """
}

process CollateData {
    input:
    path processed

    output:
    path "${nxtse_path}/filteredIntrons.Rds", emit: filtered_introns

    script:
    """
    scripts/collateData.R
    """
}

process CoveragePlots {
    input:
    path filtered_introns

    output:
    path "${ref_path}/plots_done.txt", emit: plots_done

    script:
    """
    scripts/coveragePlots.R
    touch ${ref_path}/plots_done.txt
    """
}

// Define workflow
workflow {
    // Step 1: Get resources
    getResources(config.genome_fasta, config.genome_gtf)

    // Step 2: Build references
    SpliceRef(getResources.out.copied_genome, getResources.out.copied_gtf)
    StarRef(getResources.out.copied_genome, getResources.out.copied_gtf)

    // Step 3: Align reads
    Channel.from(config.samples)
        .map { sample -> 
            tuple(
                sample, 
                file("${config.fastq_folder}/${sample}_1.fq.gz"), 
                file("${config.fastq_folder}/${sample}_2.fq.gz")
            )
        }
        .set { fastq_channel }

    StarAlignment(StarRef.out.star_done, fastq_channel)

    // Step 4: Process BAMs
    ProcessBams(StarAlignment.out.aligned_bam, SpliceRef.out.splicewiz_ref)

    // Step 5: Collate data
    CollateData(ProcessBams.out.processed.collect())

    // Step 6: Generate coverage plots
    CoveragePlots(CollateData.out.filtered_introns)

    // Final output
    CoveragePlots.out.plots_done.view { "Pipeline completed: ${it}" }
}