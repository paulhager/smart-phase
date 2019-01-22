#!/bin/bash

# STARTING CONFIGURATION
# ./CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf.gz
# ./YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf.gz
# ./shapeit.v2.904.2.6.32-696.18.7.el6.x86_64
# ./whatshap-comparison-experiments
# ./reference/
#       1000GP_Phase3.sample
#       genetic_map_chr1_combined_b37.txt
#       1000GP_Phase3_chr1.legend.gz
#       1000GP_Phase3_chr1.hap.gz
#       human_g1k_v37.fasta
capture_bed=../BED/AGV6UTR_covered_merged.bed
gene_bed=../BED/allGeneRegionsCanonical.HG19.GRCh37.bed

# DOWNLOAD LOCATIONS
ftp_vcf=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/
url_res=http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/

chromosomes=(1 19)


sorted_gene_bed=$( basename $gene_bed )
sorted_gene_bed=./reference/${sorted_gene_bed/.bed/.sort.bed}
sort -k1,1 -k2,2n $gene_bed > $sorted_gene_bed

merged_gene_bed=${sorted_gene_bed/.bed/.merge.bed}
bedtools merge -i $sorted_gene_bed > $merged_gene_bed

intersect_bed=reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.bed
bedtools intersect -a $capture_bed -b $merged_gene_bed > $intersect_bed

intersect_cut_bed=${intersect_bed/.bed/.cut.bed}
cut -c4- $intersect_bed > $intersect_cut_bed

