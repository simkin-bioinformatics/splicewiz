gtf_file = '/home/charlie/big_data/genomes_and_TE/BDGP6.32/dm_bdgp632_with_copia.gtf'
output_bed_file = '/home/charlie/big_data/genomes_and_TE/BDGP6.32/dm_bdgp632_copia_regions.bed'

with open(gtf_file, 'r') as gtf, open(output_bed_file, 'w') as bed:
    for line in gtf:
        if 'COPIA_DM_I\tgene' in line:
            line = line.strip().split()
            chrom, start, end = line[0], line[3], line[4]
            bed.write(f'{chrom}\t{start}\t{end}\n')