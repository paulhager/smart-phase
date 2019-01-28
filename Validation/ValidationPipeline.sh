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

# STARTING CONFIGURATION (non-existing files will be downloaded)
# ./CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf.gz
# ./YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf.gz
# ./reference/
#	1000GP_Phase3.sample
#	genetic_map_chr1_combined_b37.txt
#	1000GP_Phase3_chr1.legend.gz
#	1000GP_Phase3_chr1.hap.gz
#  	human_g1k_v37.fasta
# ./sim
# ./whatshap-comparison-experiments
shapeit=./bin/shapeit
art_illumina=./bin/art_illumina
capture_bed=../BED/AGV6UTR_covered_merged.bed
gene_bed=../BED/allGeneRegionsCanonical.HG19.GRCh37.bed

# DOWNLOAD LOCATIONS
ftp_vcf=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/
url_res=http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/

chromosomes=(1 19)

# check whether SHAPEIT is available
if [[ ! -f $shapeit ]]; then
	echo "SHAPEIT binary is missing! Please download it as explained in the README."
	exit
fi

# check whether ART is available
if [[ ! -f $art_illumina ]]; then
        echo "ART binary is missing! Please download it as explained in the README."
        exit
fi

# Preliminary work that only needs to be done once
#rm -r ./sim
mkdir sim
#rm -r ./reference
mkdir ./reference
#rm -r ./results
mkdir ./results
#rm -r ./whatshap-comparison-experiments
mkdir ./whatshap-comparison-experiments
mkdir ./whatshap-comparison-experiments/scripts

# download and prepare reference genome
curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz -o ./reference/human_g1k_v37.fasta.gz
gunzip ./reference/human_g1k_v37.fasta.gz
bwa index reference/human_g1k_v37.fasta
samtools faidx reference/human_g1k_v37.fasta

# download required scripts provided by WhatsHap
curl https://bitbucket.org/whatshap/phasing-comparison-experiments/raw/08cd648ea5a4d19d8efa61f9be658c914a964f3b/scripts/artificial-child.py -o ./whatshap-comparison-experiments/scripts/artificial-child.py
chmod 777 ./whatshap-comparison-experiments/scripts/artificial-child.py
curl https://bitbucket.org/whatshap/phasing-comparison-experiments/raw/08cd648ea5a4d19d8efa61f9be658c914a964f3b/scripts/genomesimulator.py -o ./whatshap-comparison-experiments/scripts/genomesimulator.py
chmod 777 ./whatshap-comparison-experiments/scripts/genomesimulator.py

gp_samples=./reference/1000GP_Phase3.sample
curl $url_res/1000GP_Phase3.sample -o $gp_samples
for chrNum in ${chromosomes[@]}; do
  curl $url_res/1000GP_Phase3_chr$chrNum.hap.gz -o ./reference/1000GP_Phase3_chr$chrNum.hap.gz
  curl $url_res/1000GP_Phase3_chr$chrNum.legend.gz -o ./reference/1000GP_Phase3_chr$chrNum.legend.gz
  curl $url_res/genetic_map_chr${chrNum}_combined_b37.txt -o ./reference/genetic_map_chr${chrNum}_combined_b37.txt
  awk 'NR==1 {{print}} NR>1 {{print $1, $2*1.5, $3*1.5 }}' reference/genetic_map_chr${chrNum}_combined_b37.txt > reference/genetic_map_chr${chrNum}.txt
done

# preparing BED file representing captured protein-coding exons (intersection of captured regions and protein-coding genes)
sorted_gene_bed=$( basename $gene_bed )
sorted_gene_bed=./reference/${sorted_gene_bed/.bed/.sort.bed}
sort -k1,1 -k2,2n $gene_bed > $sorted_gene_bed

merged_gene_bed=${sorted_gene_bed/.bed/.merge.bed}
bedtools merge -i $sorted_gene_bed > $merged_gene_bed

intersect_bed=reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.bed
bedtools intersect -a $capture_bed -b $merged_gene_bed > $intersect_bed

intersect_cut_bed=${intersect_bed/.bed/.cut.bed}
cut -c4- $intersect_bed > $intersect_cut_bed

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
    allVCFPath=$ftp_vcf/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf.gz
    ;;
  "2")
    sample=YRI
    father=NA19239
    mother=NA19238
    child=NA19240
    parents=($mother $father)
    family=($child $mother $father)
    allVCFPath=$ftp_vcf/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.vcf.gz
    ;;
  *)
    ;;
  esac

  # Clean and create structure
  rm -r ./$sample
  mkdir ./$sample

  mkdir ./sim/$sample
  mkdir -p sim/tmp

  # Download required files
  allVCFFile=${allVCFPath/*\//}
  curl $allVCFPath -o ./$sample/$allVCFFile

  # Unzip unphased vcfs
  allVCFFileUnzip=${allVCFFile/.gz/}
  gunzip $sample/$allVCFFile 

  for chrNum in ${chromosomes[@]}; do
    allVCFFileUnzipChr=${allVCFFileUnzip/.vcf/.chr$chrNum.vcf}
    # Extract chr1 and chr19 from vcf
    ""awk '/^#/ || ($1 == var)' var="$chrNum" $sample/$allVCFFileUnzip > $sample/$allVCFFileUnzipChr""

    # Zip and index
    bgzip $sample/$allVCFFileUnzipChr
    allVCFFileChr=$allVCFFileUnzipChr.gz
    tabix -p vcf $sample/$allVCFFileChr

    # Extract biallelic sites
    allVCFFileChrBiallelic=${allVCFFileChr/.vcf.gz/_biallelic.vcf.gz}
    vcftools --gzvcf $sample/$allVCFFileChr --min-alleles 2 --max-alleles 2 --recode --stdout | bgzip -c > $sample/$allVCFFileChrBiallelic
    tabix -p vcf $sample/$allVCFFileChrBiallelic

    # Run SHAPEIT check
    genmap=./reference/genetic_map_chr${chrNum}_combined_b37.txt
    refhaps=./reference/1000GP_Phase3_chr${chrNum}.hap.gz
    legend=./reference/1000GP_Phase3_chr${chrNum}.legend.gz
    $shapeit -check -V $sample/$allVCFFileChrBiallelic -M $genmap --input-ref $refhaps $legend $gp_samples --output-log ./$sample/shapeit_check_chr${chrNum} > ./$sample/shapeit_check_chr${chrNum}.log 2>&1

    # Run SHAPEIT phasing
    exclude=./$sample/shapeit_check_chr${chrNum}.snp.strand.exclude
    $shapeit -V $sample/$allVCFFileChrBiallelic --exclude-snp $exclude -M $genmap --input-ref $refhaps $legend $gp_samples -O ./$sample/shapeit_phase_chr${chrNum} > $sample/shapeit_phase_chr${chrNum}.log 2>&1

    # Convert SHAPEIT result to VCF
    allVCFPhased=${allVCFFileChrBiallelic/.vcf.gz/_phased.vcf}
    $shapeit -convert --input-haps ./$sample/shapeit_phase_chr${chrNum} --output-vcf ./$sample/$allVCFPhased > ./$sample/shapeit_convert_chr${chrNum}.log 2>&1

    allVCFPhasedZip=${allVCFPhased/.vcf/.vcf.gz}
    bgzip $sample/$allVCFPhased
    tabix -p vcf $sample/$allVCFPhasedZip

    # Split VCF into three
    for familyMember in ${family[@]}; do
      singleVCFPhased=${allVCFPhased/genotypes./genotypes.${familyMember}.}
      bcftools view -s $familyMember $sample/$allVCFPhasedZip > $sample/$singleVCFPhased
    done
   
    # Create artifical child
    motherVCF=${allVCFPhased/genotypes./genotypes.${mother}.}                                          
    fatherVCF=${allVCFPhased/genotypes./genotypes.${father}.}
    childVCF=sim.$sample.$child.chr${chrNum}.phased.vcf
    whatshap-comparison-experiments/scripts/artificial-child.py reference/genetic_map_chr${chrNum}.txt $sample/$motherVCF $sample/$fatherVCF $child sim/$sample/$sample.$father.chr${chrNum}.true.recomb sim/$sample/$sample.$mother.chr${chrNum}.true.recomb > sim/$sample/$childVCF
  
    # Create child true haplotype fastas
    whatshap-comparison-experiments/scripts/genomesimulator.py -c $chrNum sim/$sample/$childVCF reference/human_g1k_v37.fasta sim/tmp/
 
    motherVCFZip=${motherVCF/.vcf/.vcf.gz}
    bgzip $sample/$motherVCF
    tabix -p vcf $sample/$motherVCFZip
  
    fatherVCFZip=${fatherVCF/.vcf/.vcf.gz}
    bgzip $sample/$fatherVCF
    tabix -p vcf $sample/$fatherVCFZip
  
    childVCFZip=${childVCF/.vcf/.vcf.gz}
    bgzip sim/$sample/$childVCF
    tabix -p vcf sim/$sample/$childVCFZip
  
    # Merge artifical trio
    trioVCF=sim.${sample}.trio.chr${chrNum}.phased.vcf
    vcf-merge $sample/$motherVCFZip $sample/$fatherVCFZip sim/$sample/$childVCFZip > sim/$sample/$trioVCF
  
  done

  rm $sample/$allVCFFileUnzip

  # Create ped
  echo -e "${sample}\t${child}\t${father}\t${mother}\t0\t2" > reference/$sample.ped

  # Cut desired regions out of fasta and move to sample folders
  for f in sim/tmp/${child}*.fasta; do
    samtools faidx $f
    outFile=${f/.fasta/_kit.fasta}
    outFile=${outFile/tmp/$sample}
    bedtools getfasta -fi $f -bed $capture_bed -fo $outFile
    rm $f
    rm $f.fai
  done

  # Simulate reads of both haplotypes using ART
  length=150
  cv=100
  frag=400
  sdev=100

  for f in sim/$sample/*_kit.fasta; do
    sample_chr=$( basename $f | cut -d "_" -f1 )
    $art_illumina -ss HSXt -i $f -p -l $length -f $cv -m $frag -s $sdev -o sim/$sample/simulated.art.hsxt.${length}l.${cv}fc.${frag}m.${sdev}s.${sample_chr}. > sim/$sample/art_simulation_${sample_chr}.log 2>&1
  done

  # Rename headers to prevent duplicate names
  for f in sim/$sample/*.1.[12].fq; do
    out=${f/.fq/.fixed_header.fq.gz}
    cat $f | awk '{if(NR%4==1) $0="@hapl1_"$0; print;}' | gzip > $out 
  done

  for f in sim/$sample/*.2.[12].fq; do
    out=${f/.fq/.fixed_header.fq.gz}
    cat $f | awk '{if(NR%4==1) $0="@hapl2_"$0; print;}' | gzip > $out 
  done

  # Merge Haplotypes
  for f in sim/$sample/*1.1.fixed_header.fq.gz; do
    f2=${f/1.1.fixed_header.fq.gz/2.1.fixed_header.fq.gz}
    zcat $f $f2 | gzip > ${f/1.1.fixed_header.fq.gz/1.fq.gz} 
  done

  for f in sim/$sample/*1.2.fixed_header.fq.gz; do
    f2=${f/1.2.fixed_header.fq.gz/2.2.fixed_header.fq.gz}
    zcat $f $f2 | gzip > ${f/1.2.fixed_header.fq.gz/2.fq.gz} 
  done

  # Run BWA
  threads=4
  for f in sim/$sample/*1.fq.gz; do
    id=$( basename $f | cut -d "." -f8 )
    f2=${f/1.fq.gz/2.fq.gz}
    out=${f/.1.fq.gz/.bam}
    bwa mem -t $threads -R "@RG\tID:${id}\tSM:${id}" reference/human_g1k_v37.fasta $f $f2 | samtools sort -@ $threads -o $out -
    samtools index $out 
  done 

  # Extract variants in kit and intersection of kit with all gene regions canonical from USCS table browser.
  for f in sim/$sample/sim.${sample}.trio.*.vcf; do
    out=${f/.vcf/.intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37}
    vcftools --vcf $f --bed $intersect_cut_bed --out $out --recode --keep-INFO-all
    out=${out}.recode.vcf
    bgzip $out
    out=${out}.gz
    tabix -p vcf $out
  done

  # Run SmartPhase
  for chrNum in ${chromosomes[@]}; do
    vcf=./sim/$sample/sim.${sample}.trio.chr${chrNum}.phased.intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz
    bam=./sim/$sample/simulated.art.hsxt.150l.100fc.400m.100s.${child}.chr${chrNum}.bam
    out_raw=$( basename $vcf )

    # read-only phasing
    out=${out_raw/.vcf.gz/.SmartPhase_read-only.tsv}
    log=${out/.tsv/.log}
    java -jar ../smartPhase.jar -g $merged_gene_bed -a $vcf -p $child -r $bam -m 60 -o results/$out -v -x > results/$log 2>&1
    perl ./extract_SmartPhase_stats.pl results $sample $child $chrNum "read-only" 10 0.1
    
    # read and trio phasing
    out=${out_raw/.vcf.gz/.SmartPhase_read-and-trio.tsv}
    log=${out/.tsv/.log}
    java -jar ../smartPhase.jar -g $merged_gene_bed -a $vcf -p $child -r $bam -m 60 -d reference/${sample}.ped -o results/$out -v -x -t > results/$log 2>&1
    perl ./extract_SmartPhase_stats.pl results $sample $child $chrNum "read-and-trio" 10 0.1
  done

  # Run WhatsHap
  for chrNum in ${chromosomes[@]}; do
    vcf=./sim/$sample/sim.${sample}.trio.chr${chrNum}.phased.intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz
    bam=./sim/$sample/simulated.art.hsxt.150l.100fc.400m.100s.${child}.chr${chrNum}.bam
    eval_pairs=./results/sim.${sample}.trio.chr${chrNum}.phased.intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.evaluation_pairs.tsv
    out_raw=$( basename $vcf )

    # unphase simluated VCF
    uv=${vcf/.vcf.gz/_unphased.vcf.gz}
    whatshap unphase $vcf | bgzip > $uv
    tabix -p vcf $uv

    # read-only phasing
    out=${out_raw/.vcf.gz/.WhatsHap_read-only.vcf.gz}
    log=${out/.vcf.gz/.log}
    whatshap phase --mapq 60 --indels --sample $child -o results/$out $uv $bam > results/$log 2>&1
    perl ./extract_WhatsHap_results.pl results/$out $eval_pairs

    # read and trio phasing
    out=${out_raw/.vcf.gz/.WhatsHap_read-and-trio.vcf.gz}
    log=${out/.vcf.gz/.log}
    whatshap phase --mapq 60 --indels --ped reference/${sample}.ped -o results/$out $uv $bam > results/$log 2>&1
    perl ./extract_WhatsHap_results.pl results/$out $eval_pairs
  done

done
