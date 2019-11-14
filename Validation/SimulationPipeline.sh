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
# ./bin/shapeit (needs to be added by the user as explained in the README)
# ./Wessim
# ./GemSim	(needs to be added by the user as explained in the README)
shapeit=./bin/shapeit
capture_bed=../BED/AGV6UTR_covered_merged.bed
gene_bed=../BED/allGeneRegionsCanonical.HG19.GRCh37.bed

# !! The probe file cannot be shared for proprietary reasons. If you want to run Wessim2 as shown below (Step 8), you have to request it from Agilent. !!
probe_file=V6_UTR_Probes_chr1_chr19.txt

# Error model file provided by GemSim
model_file=./GemSIM_v1.6/models/ill100v5_p.gzip

# Download locations
ftp_vcf=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/
url_res=http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/

chromosomes=(1 19)

# Check whether SHAPEIT is available
if [[ ! -f $shapeit ]]; then
	echo "SHAPEIT binary is missing! Please download it as explained in the README."
	exit
fi

# Preliminary work that only needs to be done once
mkdir -p ./CEU
mkdir -p ./YRI
mkdir -p ./reference
mkdir -p ./sim/CEU
mkdir -p ./sim/YRI
mkdir -p ./whatshap-comparison-experiments/scripts

##########
# Step 1 #
##########

# download and prepare reference genome
curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz -o ./reference/human_g1k_v37.fasta.gz
gunzip ./reference/human_g1k_v37.fasta.gz
bwa index reference/human_g1k_v37.fasta
samtools faidx reference/human_g1k_v37.fasta

# download reference data
gp_samples=./reference/1000GP_Phase3.sample
curl $url_res/1000GP_Phase3.sample -o $gp_samples
for chrNum in ${chromosomes[@]}; do
  curl $url_res/1000GP_Phase3_chr${chrNum}.hap.gz -o ./reference/1000GP_Phase3_chr${chrNum}.hap.gz
  curl $url_res/1000GP_Phase3_chr${chrNum}.legend.gz -o ./reference/1000GP_Phase3_chr${chrNum}.legend.gz
  curl $url_res/genetic_map_chr${chrNum}_combined_b37.txt -o ./reference/genetic_map_chr${chrNum}_combined_b37.txt
  awk 'NR==1 {{print}} NR>1 {{print $1, $2*1.5, $3*1.5 }}' reference/genetic_map_chr${chrNum}_combined_b37.txt > reference/genetic_map_chr${chrNum}.txt
done

# download required scripts provided by WhatsHap
curl https://bitbucket.org/whatshap/phasing-comparison-experiments/raw/08cd648ea5a4d19d8efa61f9be658c914a964f3b/scripts/artificial-child.py -o ./whatshap-comparison-experiments/scripts/artificial-child.py
chmod 777 ./whatshap-comparison-experiments/scripts/artificial-child.py
curl https://bitbucket.org/whatshap/phasing-comparison-experiments/raw/08cd648ea5a4d19d8efa61f9be658c914a964f3b/scripts/genomesimulator.py -o ./whatshap-comparison-experiments/scripts/genomesimulator.py
chmod 777 ./whatshap-comparison-experiments/scripts/genomesimulator.py

# download Wessim and prepare probe file
git clone https://github.com/sak042/Wessim
python ./Wessim/Wessim_ver_1.0/Prep_Probe2Fa.py $probe_file

##########
# Step 2 #
##########

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

##########
# Step 3 #
##########

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

##########
# Step 4 #
##########

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

##########
# Step 5 #
##########

    # Split SHAPEIT phased VCF into three
    for familyMember in ${family[@]}; do
      singleVCFPhased=${allVCFPhased/genotypes./genotypes.${familyMember}.}
      bcftools view -s $familyMember $sample/$allVCFPhasedZip > $sample/$singleVCFPhased
    done

    # Create artifical child
    motherVCF=${allVCFPhased/genotypes./genotypes.${mother}.}                                          
    fatherVCF=${allVCFPhased/genotypes./genotypes.${father}.}
    childVCF=sim.$sample.$child.chr${chrNum}.phased.vcf
    whatshap-comparison-experiments/scripts/artificial-child.py reference/genetic_map_chr${chrNum}.txt $sample/$motherVCF $sample/$fatherVCF $child sim/$sample/$sample.$father.chr${chrNum}.true.recomb sim/$sample/$sample.$mother.chr${chrNum}.true.recomb > sim/$sample/$childVCF

##########
# Step 6 #
##########
  
    # Create child true haplotype fastas
    whatshap-comparison-experiments/scripts/genomesimulator.py -c $chrNum sim/$sample/$childVCF reference/human_g1k_v37.fasta sim/$sample/

##########
# Step 7 #
##########
 
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

    # Extract variants in kit and intersection of kit with all gene regions canonical from USCS table browser.
    cutVCF=sim/$sample/${trioVCF/.vcf/.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37}
    vcftools --vcf sim/$sample/$trioVCF --bed $intersect_cut_bed --out $cutVCF --recode --keep-INFO-all
    cutVCF=${cutVCF}.recode.vcf
    mv $cutVCF ${cutVCF/.recode/}
    cutVCF=${cutVCF/.recode/}
    bgzip $cutVCF
    cutVCF=${cutVCF}.gz
    tabix -p vcf $cutVCF

  done

  rm $sample/$allVCFFileUnzip

##########
# Step 8 #
##########

  # Simulate reads of both haplotypes using Wessim2
  length=100
  frag=250
  sdev=100

  for f in sim/$sample/*.fasta; do
    sample_chr=$( basename $f | cut -d "." -f1-3 )
    samtools $f
    faToTwoBit $f ${f/.fasta/.2bit}
  done

  cd ./Wessim/Wessim_ver_1.0

  for f in ../../sim/$sample/*2bit; do
    sampleChrHap=$( basename $f | cut -d "." -f1-3)
    chr=$( echo $sampleChrHap | cut -d "." -f2 )
    if [ "$chr" == "chr1" ]; then
      reads=4500000
    else 
      reads=2500000
    fi
    out=../../sim/$sample/simulated.Wessim2.ill100v5_p.${length}l.${reads}reads.${frag}f.${sdev}d.${sampleChrHap}

    # Establish local blat server in background
    gfServer start -canStop localhost 6666 $( realpath $f ) &
    proc_id=$( echo $! )
    sleep 10

    # Run blat search to generate the match list
    python Prep_BlatSearch.py -R $( basename $f ) -P ${probe_file}.fa -o ${f/.2bit/.psl}

    # Run Wessim2 in probe hybridization mode
    threads=16
    python Wessim2.py -R ${f/.2bit/.fasta} -P ${probe_file}.fa -B ${f/.2bit/.psl} -f $frag -d $sdev -p -n $reads -l $length -M $model_file -t $threads -o $out -z -v > ../../sim/$sample/${sampleChrHap}_Wessim2.out 2> ../../sim/$sample/${sampleChrHap}_Wessim2.err

    # Kill local blat server
    kill $proc_id
  done

  cd ../../

  # Rename headers to prevent duplicate names
  # FASTQs of haplotype 1
  for f in sim/$sample/*.1_[12].fastq.gz; do
    out=${f/.fastq.gz/.fixed_header.fq.gz}
    cat $f | awk '{if(NR%4==1) $0="@hapl1_"$0; print;}' | gzip > $out 
  done

  # FASTQs of haplotype 2
  for f in sim/$sample/*.2_[12].fastq.gz; do
    out=${f/.fastq.gz/.fixed_header.fq.gz}
    cat $f | awk '{if(NR%4==1) $0="@hapl2_"$0; print;}' | gzip > $out 
  done

##########
# Step 9 #
##########

  # Merge Haplotypes
  # forward reads
  for f in sim/$sample/*1_1.fixed_header.fq.gz; do
    f2=${f/1_1.fixed_header.fq.gz/2_1.fixed_header.fq.gz}
    zcat $f $f2 | gzip > ${f/1_1.fixed_header.fq.gz/1.fq.gz} 
  done

  # reverse reads
  for f in sim/$sample/*1_2.fixed_header.fq.gz; do
    f2=${f/1_2.fixed_header.fq.gz/2_2.fixed_header.fq.gz}
    zcat $f $f2 | gzip > ${f/1_2.fixed_header.fq.gz/2.fq.gz} 
  done

done
