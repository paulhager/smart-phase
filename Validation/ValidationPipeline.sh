Pipeline Commands NA12878

DOWNLOAD LOCATIONS:
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/

STARTING CONFIGURATION

.//$sample
    ..//$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf.gz

.//shapeit.v2.904.2.6.32-696.18.7.el6.x86_64
.//whatshap-comparison-experiments
.//reference
	..//1000GP_Phase3.sample
	..//genetic_map_chr1_combined_b37.txt
	..//1000GP_Phase3_chr1.legend.gz
	..//1000GP_Phase3_chr1.hap.gz
  ..//human_g1k_v37.fasta
	..//AGV6UTR_covered_merged.bed # Don't know where Tim got this from
  ..//allGeneRegionsCanonical.HG19.GRCh37.sort.merge.bed # Can't be downloaded automatically

chromosomes=(1 19)

# Preliminary work that only needs to be done once
rm -r ./sim
mkdir sim
rm -r ./reference
mkdir ./reference
rm -r ./whatshap-comparison-experiments
mkdir ./whatshap-comparison-experiments
mkdir ./whatshap-comparison-experiments/scripts
curl http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3.sample -o ./reference/1000GP_Phase3.sample
curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz -o ./reference/human_g1k_v37.fasta.gz
gunzip -c ./reference/human_g1k_v37.fasta.gz > ./reference/human_g1k_v37.fasta
samtools faidx reference/human_g1k_v37.fasta
curl https://bitbucket.org/whatshap/phasing-comparison-experiments/raw/08cd648ea5a4d19d8efa61f9be658c914a964f3b/scripts/artificial-child.py -o ./whatshap-comparison-experiments/scripts/artificial-child.py
chmod 777 ./whatshap-comparison-experiments/scripts/artificial-child.py
curl https://bitbucket.org/whatshap/phasing-comparison-experiments/raw/08cd648ea5a4d19d8efa61f9be658c914a964f3b/scripts/genomesimulator.py -o ./whatshap-comparison-experiments/scripts/genomesimulator.py
chmod 777 ./whatshap-comparison-experiments/scripts/genomesimulator.py
for chrNum in ${chromosomes[@]}; do
  curl http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3_chr$chrNum.hap.gz -o ./reference/1000GP_Phase3_chr$chrNum.hap.gz
  curl http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3_chr$chrNum.legend.gz -o ./reference/1000GP_Phase3_chr$chrNum.legend.gz
  curl http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/genetic_map_chr${chrNum}_combined_b37.txt -o ./reference/genetic_map_chr${chrNum}_combined_b37.txt
  awk 'NR==1 {{print}} NR>1 {{print $1, $2*1.5, $3*1.5 }}' reference/genetic_map_chr${chrNum}_combined_b37.txt > reference/genetic_map_chr${chrNum}.txt
done

bedtools intersect -a reference/AGV6UTR_covered_merged.bed -b reference/allGeneRegionsCanonical.HG19.GRCh37.sort.merge.bed > reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.bed

cut -c4- reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.bed > reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.cut.bed

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
    allVCFPath=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf.gz
    ;;
  "2")
    sample=YRI
    parents=(NA19238 NA19239)
    father=NA19239
    mother=NA19238
    child=NA19240
    parents=($mother $father)
    family=($child $mother $father)
    allVCFPath=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.vcf.gz
    ;;
  *)
    ;;
  esac

  # Clean and create structure
  rm -r ./$sample
  mkdir ./$sample

  mkdir ./sim/$sample

  # Download required files
  allVCFFile=${allVCFPath/*\//}
  curl $allVCFPath -o ./$sample/$allVCFFile

  # Unzip unphased vcfs
  allVCFFileUnzip=${allVCFFile/.gz/}
  gunzip -c $sample/$allVCFFile > $sample/$allVCFFileUnzip 

  for chrNum in ${chromosomes[@]}; do
    # Extract chr1 and chr19 from vcf
    ""awk '/^#/ || ($1 == var)' var="$chrNum" $sample/$sample.wgs.consensus.20131118.snps_indels.*.genotypes.vcf > $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr$chrNum.vcf""

    # Zip and index
    bgzip -c $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr$chrNum.vcf > $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr$chrNum.vcf.gz
    tabix $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr$chrNum.vcf.gz
  done

#@@@@@@ TIM NEEDS TO DO THIS PART BECAUSE SHAPEIT DOES NOT WORK ON OSX
  # Run SHAPEIT Check, then SHAPEIT, then convert to vcf
  shapeit -check -V $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1.vcf -M reference/genetic_map_chr1_combined_b37.txt --input-ref reference/1000GP_Phase3_chr1.hap.gz reference/1000GP_Phase3_chr1.legend.gz reference/1000GP_Phase3.sample --output-log shapeit/NA12878.chr1 || true) > shapeit/NA12878.trio.chr1.check.log 2>&1
  shapeit -V $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1.vcf --exclude-snp shapeit/NA12878.trio.chr1.snp.strand.exclude' -M reference/genetic_map_chr1_combined_b37.txt --input-ref reference/1000GP_Phase3_chr1.hap.gz reference/1000GP_Phase3_chr1.legend.gz reference/1000GP_Phase3.sample -O $sample/NA12878.trio.chr1 > shapeit/NA12878.trio.chr1.run.log 2>&
  shapeit -convert --input-haps shapeit/NA12878.trio.chr1 --output-vcf $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1_biallelic_phased.vcf > shapeit/NA12878.trio.chr1.phased.vcf.log 2>&1

  shapeit -check -V $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19.vcf -M reference/genetic_map_chr19_combined_b37.txt --input-ref reference/1000GP_Phase3_chr19.hap.gz reference/1000GP_Phase3_chr19.legend.gz reference/1000GP_Phase3.sample --output-log shapeit/NA12878.chr19 || true) > shapeit/NA12878.trio.chr19.check.log 2>&1
  shapeit -V $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19.vcf --exclude-snp shapeit/NA12878.trio.chr19.snp.strand.exclude' -M reference/genetic_map_chr19_combined_b37.txt --input-ref reference/1000GP_Phase3_chr19.hap.gz reference/1000GP_Phase3_chr19.legend.gz reference/1000GP_Phase3.sample -O $sample/NA12878.trio.chr19 > shapeit/NA12878.trio.chr19.run.log 2>&
  shapeit -convert --input-haps shapeit/NA12878.trio.chr19 --output-vcf $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19_biallelic_phased.vcf > shapeit/NA12878.trio.chr19.phased.vcf.log 2>&1
#@@@@@@ TIM NEEDS TO DO THIS PART BECAUSE SHAPEIT DOES NOT WORK ON OSX

  # Split VCF into three
  for chrNum in ${chromosomes[@]}; do
    for familyMember in ${family[@]}; do
      bcftools view -s $familyMember $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr${chrNum}_biallelic_phased.vcf.gz > $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.$familyMember.chr${chrNum}_biallelic_phased.vcf
    done
  done

  # Create artifical child
  for chrNum in ${chromosomes[@]}; do
    whatshap-comparison-experiments/scripts/artificial-child.py reference/genetic_map_chr${chrNum}.txt $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.$mother.chr${chrNum}_biallelic_phased.vcf $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.$father.chr${chrNum}_biallelic_phased.vcf $child sim/$sample/$sample.$father.chr${chrNum}.true.recomb sim/$sample/$sample.$mother.chr${chrNum}.true.recomb > sim/$sample/sim.$sample.$child.chr${chrNum}.phased.vcf
  done

  # Merge artifical trio
  for chrNum in ${chromosomes[@]}; do
    for parent in ${parents[@]}; do
      bgzip -c $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.$parent.chr${chrNum}_biallelic_phased.vcf > $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.$parent.chr${chrNum}_biallelic_phased.vcf.gz
      tabix $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.$parent.chr${chrNum}_biallelic_phased.vcf.gz
    done

    bgzip -c sim/$sample/sim.$sample.$child.chr${chrNum}.phased.vcf > sim/$sample/sim.$sample.$child.chr${chrNum}.phased.vcf.gz
    tabix sim/$sample/sim.$sample.$child.chr${chrNum}.phased.vcf.gz

    vcf-merge $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.$mother.chr${chrNum}_biallelic_phased.vcf.gz $sample/$sample.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.$father.chr${chrNum}_biallelic_phased.vcf.gz sim/$sample/sim.$sample.$child.chr${chrNum}.phased.vcf.gz > sim/$sample/sim.$sample.trio.chr${chrNum}.phased.vcf
  done

  # Create child true haplotype fastas
  mkdir -p sim/tmp

  for chrNum in ${chromosomes[@]}; do
    whatshap-comparison-experiments/scripts/genomesimulator.py -c $chrNum sim/$sample/sim.$sample.$child.chr$chrNum.phased.vcf reference/human_g1k_v37.fasta sim/tmp/
  done

  # Create ped
  echo $sample	$child	$father	$mother	0	2 > reference/$sample.ped

  # Cut desired regions out of fasta
  for f in sim/tmp/*fasta; do
    outFile=${f/.fasta/_kit.fasta}
    outFile=${outFile/tmp/$sample}
    bedtools getfasta -fi $f -bed reference/AGV6UTR_covered_merged.bed -fo $outFile
    rm $f
    rm $f.fai
  done

  # Simulate reads of both haplotypes using ART
  art=art_bin_MountRainier/art_illumina
  length=150
  cv=100
  frag=400
  sdev=100

  for f in sim/$sample/*_kit.fasta; do
    sample=$( basename $f | cut -d "_" -f1 )
    echo start for $sample
    $art -ss HSXt -i $f -p -l $length -f $cv -m $frag -s $sdev -o sim/$sample/simulated.art.hsxt.${length}l.${cv}fc.${frag}m.${sdev}s.${sample}. > ${sample}.out 2> ${sample}.err &
  done

  # Rename headers to prevent duplicate names
  for f in sim/$sample/*.1.[12].fq; do
    out=${f/.fq/.fixed_header.fq.gz}
    echo $f to $out
    cat $f | awk '{if(NR%4==1) $0="@hapl1_"$0; print;}' | gzip > $out &
  done

  for f in sim/$sample/*.2.[12].fq; do
    out=${f/.fq/.fixed_header.fq.gz}
    echo $f to $out
    cat $f | awk '{if(NR%4==1) $0="@hapl2_"$0; print;}' | gzip > $out &
  done

  # Merge Haplotypes
  for f in sim/$sample/*1.1.fixed_header.fq.gz; do
    f2=${f/1.1.fixed_header.fq.gz/2.1.fixed_header.fq.gz}
    gzcat $f $f2 | gzip > ${f/1.1.fixed_header.fq.gz/1.fq.gz} &
  done

  for f in sim/$sample/*1.2.fixed_header.fq.gz; do
    f2=${f/1.2.fixed_header.fq.gz/2.2.fixed_header.fq.gz}
    gzcat $f $f2 | gzip > ${f/1.2.fixed_header.fq.gz/2.fq.gz} &
  done

  # Map using BWA
  bwa index reference/human_g1k_v37.fasta

  threads=4
  for f in sim/$sample/*1.fq.gz; do
    echo $f
    id=$( basename $f | cut -d "." -f8 )
    echo $id
    f2=${f/1.fq.gz/2.fq.gz}
    echo $f2
    out=${f/.1.fq.gz/.bam}
    echo $out
    bwa mem -t $threads -R "@RG\tID:${id}\tSM:${id}" reference/human_g1k_v37.fasta $f $f2 | samtools sort -@ $threads -o $out -
    samtools index $out &
  done 


  # Extract variants in kit and intersection of kit with all gene regions canonical from USCS table browser.
  for f in sim/$sample/sim.$sample.trio.high_coverage_pcr_free.phased.vcf; do
    out=${f/.vcf/.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37}

    vcftools --vcf $f --bed reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.cut.bed --out $out --recode --keep-INFO-all

    out=$out.recode.vcf

    bgzip -c $out > $out.gz

    tabix $out.gz
  done

  # Run SmartPhase
  java -jar ~/smart-phase/smartPhase.jar -g reference/allGeneRegionsCanonical.HG19.GRCh37.bed -a sim/$sample/sim.$sample.trio.chr1.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz -p NA12878 -r sim/$sample/simulated.art.hsxt.150l.100fc.400m.100s.NA12878.chr1.bam -m 60 -d reference/$sample.ped -o results/smartPhase.sim.NA12878.chr1.NOtrio.AGV6UTR.allGeneRegionsCanonical.results.tsv -v -x

  java -jar ~/smart-phase/smartPhase.jar -g reference/allGeneRegionsCanonical.HG19.GRCh37.bed -a sim/$sample/sim.$sample.trio.chr1.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz -p NA12878 -r sim/$sample/simulated.art.hsxt.150l.100fc.400m.100s.NA12878.chr1.bam -m 60 -d reference/$sample.ped -o results/smartPhase.sim.NA12878.chr1.trio.AGV6UTR.allGeneRegionsCanonical.results.tsv -v -x -t

  java -jar ~/smart-phase/smartPhase.jar -g reference/allGeneRegionsCanonical.HG19.GRCh37.bed -a sim/$sample/sim.$sample.trio.chr19.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz -p NA12878 -r sim/$sample/simulated.art.hsxt.150l.100fc.400m.100s.NA12878.chr19.bam -m 60 -d reference/$sample.ped -o results/smartPhase.sim.NA12878.chr19.NOtrio.AGV6UTR.allGeneRegionsCanonical.results.tsv -v -x

  java -jar ~/smart-phase/smartPhase.jar -g reference/allGeneRegionsCanonical.HG19.GRCh37.bed -a sim/$sample/sim.$sample.trio.chr19.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz -p NA12878 -r sim/$sample/simulated.art.hsxt.150l.100fc.400m.100s.NA12878.chr19.bam -m 60 -d reference/$sample.ped -o results/smartPhase.sim.NA12878.chr19.trio.AGV6UTR.allGeneRegionsCanonical.results.tsv -v -x -t

done