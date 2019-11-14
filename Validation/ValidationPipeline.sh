#!/bin/bash

#  Copyright (C) 2018 the SmartPhase contributors.
#  Website: https://github.com/paulhager/smart-phase
#
#  This file is part of the SmartPhase phasing tool.
#
#  The SmartPhase phasing tool is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Download reference genome if not yet available
if [[ ! -f ./reference/human_g1k_v37.fasta.gz ]]; then
	mkdir ./reference 
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz -o ./reference/human_g1k_v37.fasta.gz
	gunzip ./reference/human_g1k_v37.fasta.gz
	bwa index reference/human_g1k_v37.fasta
	samtools faidx reference/human_g1k_v37.fasta
fi

# Prepare BED file
gene_bed=../BED/allGeneRegionsCanonical.HG19.GRCh37.bed
sorted_gene_bed=$( basename $gene_bed )
sorted_gene_bed=./reference/${sorted_gene_bed/.bed/.sort.bed}
sort -k1,1 -k2,2n $gene_bed > $sorted_gene_bed

merged_gene_bed=${sorted_gene_bed/.bed/.merge.bed}
bedtools merge -i $sorted_gene_bed > $merged_gene_bed


mkdir -p ./results
chromosomes=(1 19)

sp_outtable=./results/SP_overview.tsv
echo -e "Family\tTrio Mode?\tChr\tRuntime SP\tGenes\tVariants\tPairs\tInnocuous SP\tCis SP\tTrans SP\tPhased SP\tCleared SP\tNot phased SP\tNot found SP\tLow confidence SP\tErrors SP\tOverlap Low confidence and Errors SP" > $sp_outtable

wh_outtable=./results/WH_overview.tsv
echo -e "Family\tTrio Mode?\tChr\tRuntime WH\tCis WH\tTrans WH\tPhased WH\tNot phased WH\tErrors WH" > $wh_outtable

for iteration in 1 2; do
  # Initialize global variables
  case "$iteration" in
  "1")
    sample=CEU
    father=NA12891
    mother=NA12892
    child=NA12878
    parents=($mother $father)
    family=($child $mother $father)
    ;;
  "2")
    sample=YRI
    father=NA19239
    mother=NA19238
    child=NA19240
    parents=($mother $father)
    family=($child $mother $father)
    ;;
  *)
    ;;
  esac

  # Run BWA
  threads=16
  for f in sim/$sample/*1.fq.gz; do
    f2=${f/1.fq.gz/2.fq.gz}
    out=${f/.1.fq.gz/.bam}
    bwa mem -t $threads -R "@RG\tID:${child}\tSM:${child}" reference/human_g1k_v37.fasta $f $f2 | samtools sort -@ $threads -o $out -
    samtools index $out 
  done 

  # Run SmartPhase
  for chrNum in ${chromosomes[@]}; do
    vcf=./sim/$sample/sim.${sample}.trio.chr${chrNum}.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.vcf.gz
    bam=./sim/$sample/simulated.Wessim2.ill100v5_p.100l.2500000reads.250f.100d.${child}.chr${chrNum}.bam
    if [ "$chrNum" == "1" ]; then
          bam=./sim/$sample/simulated.Wessim2.ill100v5_p.100l.4500000reads.250f.100d.${child}.chr${chrNum}.bam
    fi 
    out_raw=$( basename $vcf )

    # read-only phasing
    out=${out_raw/.vcf.gz/.SmartPhase_read-only.tsv}
    log=${out/.tsv/.log}
    java -jar ../smartPhase.jar -g $merged_gene_bed -a $vcf -p $child -r $bam -m 60 -o results/$out -v -x > results/$log 2>&1
    perl ./extract_SmartPhase_stats.pl ./results $sample $child $chrNum "read-only" 10 0.34 >> $sp_outtable
    
    # read and trio phasing
    out=${out_raw/.vcf.gz/.SmartPhase_read+trio.tsv}
    log=${out/.tsv/.log}
    java -jar ../smartPhase.jar -g $merged_gene_bed -a $vcf -p $child -r $bam -m 60 -d ./ped/${sample}.ped -o results/$out -v -x -t > results/$log 2>&1
    perl ./extract_SmartPhase_stats.pl ./results $sample $child $chrNum "read+trio" 10 0.34 >> $sp_outtable
    
  done

  # Run WhatsHap
  for chrNum in ${chromosomes[@]}; do
    vcf=./sim/$sample/sim.${sample}.trio.chr${chrNum}.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.vcf.gz
    bam=./sim/$sample/simulated.Wessim2.ill100v5_p.100l.4500000reads.250f.100d.${child}.chr${chrNum}.bam
    eval_pairs=./results/sim.${sample}.trio.chr${chrNum}.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37_evaluation_pairs.tsv
    out_raw=$( basename $vcf )

    # unphase simluated VCF
    uv=${vcf/.vcf.gz/_unphased.vcf.gz}
    whatshap unphase $vcf | bgzip > $uv
    tabix -p vcf $uv

    # read-only phasing
    out=${out_raw/.vcf.gz/_WhatsHap_read-only.vcf.gz}
    log=./results/${out/.vcf.gz/.log}
    whatshap phase --mapq 60 --indels --sample $child -o results/$out $uv $bam > $log 2>&1
    perl ./extract_WhatsHap_results.pl results/$out $eval_pairs
    vartsv=./results/${out/.vcf.gz/_phased.tsv}
    wherr=./results/${out/.vcf.gz/_errCalls.txt}
    python ./get_phasing_errors.py $vcf results/${out/.vcf.gz/.tsv} $child > $wherr
   
    mode="no"
    cis_wh=$( grep -cP "\t1$" $vartsv )
    trans_wh=$( grep -cP "\t2$" $vartsv )
    phased_wh=$(( cis_wh + trans_wh ))
    not_phased_wh=$( grep -cP "\t4$" $vartsv )
    errors_wh=$( tail -n +2 $wherr | wc -l  | cut -d " " -f1 )
    runtime_wh=$( tail -n1 $log | rev | cut -d " " -f2 | rev )
    echo -e "${sample}\t$mode\t$chrNum\t$runtime_wh\t$cis_wh\t$trans_wh\t$phased_wh\t$not_phased_wh\t$errors_wh" >> $wh_outtable

    # read and trio phasing
    out=${out_raw/.vcf.gz/.WhatsHap_read-and-trio.vcf.gz}
    log=./results/${out/.vcf.gz/.log}
    whatshap phase --mapq 60 --indels --ped ./ped/${sample}.ped -o results/$out $uv $bam > $log 2>&1
    perl ./extract_WhatsHap_results.pl results/$out $eval_pairs
    vartsv=./results/${out/.vcf.gz/_phased.tsv}
    wherr=./results/${out/.vcf.gz/_errCalls.txt}
    python ./get_phasing_errors.py $vcf results/${out/.vcf.gz/.tsv} $child > $wherr

    mode="yes"
    cis_wh=$( grep -cP "\t1$" $vartsv )
    trans_wh=$( grep -cP "\t2$" $vartsv )
    phased_wh=$(( cis_wh + trans_wh ))
    not_phased_wh=$( grep -cP "\t4$" $vartsv )
    errors_wh=$( tail -n +2 $wherr | wc -l  | cut -d " " -f1 )
    runtime_wh=$( tail -n1 $log | rev | cut -d " " -f2 | rev )
    echo -e "${sample}\t$mode\t$chrNum\t$runtime_wh\t$cis_wh\t$trans_wh\t$phased_wh\t$not_phased_wh\t$errors_wh" >> $wh_outtable
  done

done
