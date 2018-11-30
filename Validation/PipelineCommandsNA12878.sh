Pipeline Commands NA12878

DOWNLOAD LOCATIONS:
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/

STARTING CONFIGURATION

.//NA12878
    ..//CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf.gz

.//shapeit.v2.904.2.6.32-696.18.7.el6.x86_64
.//whatshap-comparison-experiments
.//reference
	..//1000GP_Phase3.sample
	..//genetic_map_chr1_combined_b37.txt
	..//1000GP_Phase3_chr1.legend.gz
	..//1000GP_Phase3_chr1.hap.gz
	..//AGV6UTR_covered_merged.bed

# Unzip unphased vcfs
gunzip -c CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf.gz > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf

# Extract chr1 from vcf
""awk '/^#/ || ($1 == "1")' CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1.vcf""

""awk '/^#/ || ($1 == "19")' CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.vcf > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19.vcf""

# Zip and index
bgzip -c CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1.vcf > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1.vcf.gz
tabix CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1.vcf.gz

bgzip -c CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19.vcf > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19.vcf.gz
tabix CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19.vcf.gz

# Run SHAPEIT Check, then SHAPEIT, then convert to vcf
shapeit -check -V CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1.vcf -M reference/genetic_map_chr1_combined_b37.txt --input-ref reference/1000GP_Phase3_chr1.hap.gz reference/1000GP_Phase3_chr1.legend.gz reference/1000GP_Phase3.sample --output-log shapeit/NA12878.chr1 || true) > shapeit/NA12878.trio.chr1.check.log 2>&1
shapeit -V CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1.vcf --exclude-snp shapeit/NA12878.trio.chr1.snp.strand.exclude' -M reference/genetic_map_chr1_combined_b37.txt --input-ref reference/1000GP_Phase3_chr1.hap.gz reference/1000GP_Phase3_chr1.legend.gz reference/1000GP_Phase3.sample -O CEU/NA12878.trio.chr1 > shapeit/NA12878.trio.chr1.run.log 2>&
shapeit -convert --input-haps shapeit/NA12878.trio.chr1 --output-vcf CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1_biallelic_phased.vcf > shapeit/NA12878.trio.chr1.phased.vcf.log 2>&1

shapeit -check -V CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19.vcf -M reference/genetic_map_chr19_combined_b37.txt --input-ref reference/1000GP_Phase3_chr19.hap.gz reference/1000GP_Phase3_chr19.legend.gz reference/1000GP_Phase3.sample --output-log shapeit/NA12878.chr19 || true) > shapeit/NA12878.trio.chr19.check.log 2>&1
shapeit -V CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19.vcf --exclude-snp shapeit/NA12878.trio.chr19.snp.strand.exclude' -M reference/genetic_map_chr19_combined_b37.txt --input-ref reference/1000GP_Phase3_chr19.hap.gz reference/1000GP_Phase3_chr19.legend.gz reference/1000GP_Phase3.sample -O CEU/NA12878.trio.chr19 > shapeit/NA12878.trio.chr19.run.log 2>&
shapeit -convert --input-haps shapeit/NA12878.trio.chr19 --output-vcf CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19_biallelic_phased.vcf > shapeit/NA12878.trio.chr19.phased.vcf.log 2>&1

# Split VCF into three
bcftools view -s NA12878 CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1_biallelic_phased.vcf.gz > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12878.chr1_biallelic_phased.vcf
bcftools view -s NA12891 CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1_biallelic_phased.vcf.gz > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr1_biallelic_phased.vcf
bcftools view -s NA12892 CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr1_biallelic_phased.vcf.gz > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr1_biallelic_phased.vcf

bcftools view -s NA12878 CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19_biallelic_phased.vcf.gz > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12878.chr19_biallelic_phased.vcf
bcftools view -s NA12891 CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19_biallelic_phased.vcf.gz > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr19_biallelic_phased.vcf
bcftools view -s NA12892 CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.chr19_biallelic_phased.vcf.gz > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr19_biallelic_phased.vcf

# Multiply distances by factor (What is the point of this?)
awk 'NR==1 {{print}} NR>1 {{print $1, $2*1.5, $3*1.5 }}' reference/genetic_map_chr1_combined_b37.txt > reference/genetic_map_chr1.txt

awk 'NR==1 {{print}} NR>1 {{print $1, $2*1.5, $3*1.5 }}' reference/genetic_map_chr19_combined_b37.txt > reference/genetic_map_chr19.txt

# Create artifical child
whatshap-comparison-experiments/scripts/artificial-child.py reference/genetic_map_chr1.txt CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr1_biallelic_phased.vcf CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr1_biallelic_phased.vcf NA12878 sim/CEU/CEU.NA12892.chr1.true.recomb sim/CEU/CEU.NA12891.chr1.true.recomb > sim/CEU/sim.CEU.NA12878.chr1.phased.vcf

whatshap-comparison-experiments/scripts/artificial-child.py reference/genetic_map_chr19.txt CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr19_biallelic_phased.vcf CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr19_biallelic_phased.vcf NA12878 sim/CEU/CEU.NA12892.chr19.true.recomb sim/CEU/CEU.NA12891.chr19.true.recomb > sim/CEU/sim.CEU.NA12878.chr19.phased.vcf

# Merge artifical trio
bgzip -c CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr1_biallelic_phased.vcf > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr1_biallelic_phased.vcf.gz
tabix CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr1_biallelic_phased.vcf.gz
bgzip -c CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr1_biallelic_phased.vcf > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr1_biallelic_phased.vcf.gz
tabix CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr1_biallelic_phased.vcf.gz
bgzip -c sim/CEU/sim.CEU.NA12878.chr1.phased.vcf > sim/CEU/sim.CEU.NA12878.chr1.phased.vcf.gz
tabix sim/CEU/sim.CEU.NA12878.chr1.phased.vcf.gz

vcf-merge CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr1_biallelic_phased.vcf.gz CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr1_biallelic_phased.vcf.gz sim/CEU/sim.CEU.NA12878.chr1.phased.vcf.gz > sim/CEU/sim.CEU.trio.chr1.phased.vcf


bgzip -c CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr19_biallelic_phased.vcf > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr19_biallelic_phased.vcf.gz
tabix CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr19_biallelic_phased.vcf.gz
bgzip -c CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr19_biallelic_phased.vcf > CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr19_biallelic_phased.vcf.gz
tabix CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr19_biallelic_phased.vcf.gz
bgzip -c sim/CEU/sim.CEU.NA12878.chr19.phased.vcf > sim/CEU/sim.CEU.NA12878.chr19.phased.vcf.gz
tabix sim/CEU/sim.CEU.NA12878.chr19.phased.vcf.gz

vcf-merge CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12892.chr19_biallelic_phased.vcf.gz CEU/CEU.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free_v2.genotypes.NA12891.chr19_biallelic_phased.vcf.gz sim/CEU/sim.CEU.NA12878.chr19.phased.vcf.gz > sim/CEU/sim.CEU.trio.chr19.phased.vcf

# Create child true haplotype fastas
mkdir -p sim/tmp

whatshap-comparison-experiments/scripts/genomesimulator.py -c 1 sim/CEU/sim.CEU.NA12878.chr1.phased.vcf reference/human_g1k_v37.fasta sim/tmp/

whatshap-comparison-experiments/scripts/genomesimulator.py -c 19 sim/CEU/sim.CEU.NA12878.chr19.phased.vcf reference/human_g1k_v37.fasta sim/tmp/

# Create ped
echo CEU	NA12878	NA12891	NA12892	0	2 > reference/CEU.ped

# Cut desired regions out of fasta
for f in sim/tmp/*fasta; do
  outFile=${f/.fasta/_kit.fasta}
  outFile=${outFile/tmp/CEU}
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

for f in sim/CEU/*_kit.fasta; do
  sample=$( basename $f | cut -d "_" -f1 )
  echo start for $sample
  $art -ss HSXt -i $f -p -l $length -f $cv -m $frag -s $sdev -o sim/CEU/simulated.art.hsxt.${length}l.${cv}fc.${frag}m.${sdev}s.${sample}. > ${sample}.out 2> ${sample}.err &
done

# Rename headers to prevent duplicate names
for f in sim/CEU/*.1.[12].fq; do
  out=${f/.fq/.fixed_header.fq.gz}
  echo $f to $out
  cat $f | awk '{if(NR%4==1) $0="@hapl1_"$0; print;}' | gzip > $out &
done

for f in sim/CEU/*.2.[12].fq; do
  out=${f/.fq/.fixed_header.fq.gz}
  echo $f to $out
  cat $f | awk '{if(NR%4==1) $0="@hapl2_"$0; print;}' | gzip > $out &
done

# Merge Haplotypes
for f in sim/CEU/*1.1.fixed_header.fq.gz; do
  f2=${f/1.1.fixed_header.fq.gz/2.1.fixed_header.fq.gz}
  gzcat $f $f2 | gzip > ${f/1.1.fixed_header.fq.gz/1.fq.gz} &
done

for f in sim/CEU/*1.2.fixed_header.fq.gz; do
  f2=${f/1.2.fixed_header.fq.gz/2.2.fixed_header.fq.gz}
  gzcat $f $f2 | gzip > ${f/1.2.fixed_header.fq.gz/2.fq.gz} &
done

# Map using BWA
bwa index reference/human_g1k_v37.fasta

threads=4
for f in sim/CEU/*1.fq.gz; do
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


# Extract variants in kit and intersection of kit with all gene regions canonical from USCS table browser
bedtools intersect -a reference/AGV6UTR_covered_merged.bed -b reference/allGeneRegionsCanonical.HG19.GRCh37.sort.merge.bed > reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.bed

cut -c4- reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.bed > reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.cut.bed

for f in sim/CEU/sim.CEU.trio.*.phased.vcf; do
  out=${f/.vcf/.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37}

  vcftools --vcf $f --bed reference/intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.cut.bed --out $out --recode --keep-INFO-all

  out=$out.recode.vcf

  bgzip -c $out > $out.gz

  tabix $out.gz
done

# Run SmartPhase

java -jar ~/smart-phase/smartPhase.jar -g reference/allGeneRegionsCanonical.HG19.GRCh37.bed -a sim/CEU/sim.CEU.trio.chr1.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz -p NA12878 -r sim/CEU/simulated.art.hsxt.150l.100fc.400m.100s.NA12878.chr1.bam -m 60 -d reference/CEU.ped -o results/smartPhase.sim.NA12878.chr1.NOtrio.AGV6UTR.allGeneRegionsCanonical.results.tsv -v -x

java -jar ~/smart-phase/smartPhase.jar -g reference/allGeneRegionsCanonical.HG19.GRCh37.bed -a sim/CEU/sim.CEU.trio.chr1.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz -p NA12878 -r sim/CEU/simulated.art.hsxt.150l.100fc.400m.100s.NA12878.chr1.bam -m 60 -d reference/CEU.ped -o results/smartPhase.sim.NA12878.chr1.trio.AGV6UTR.allGeneRegionsCanonical.results.tsv -v -x -t

java -jar ~/smart-phase/smartPhase.jar -g reference/allGeneRegionsCanonical.HG19.GRCh37.bed -a sim/CEU/sim.CEU.trio.chr19.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz -p NA12878 -r sim/CEU/simulated.art.hsxt.150l.100fc.400m.100s.NA12878.chr19.bam -m 60 -d reference/CEU.ped -o results/smartPhase.sim.NA12878.chr19.NOtrio.AGV6UTR.allGeneRegionsCanonical.results.tsv -v -x

java -jar ~/smart-phase/smartPhase.jar -g reference/allGeneRegionsCanonical.HG19.GRCh37.bed -a sim/CEU/sim.CEU.trio.chr19.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.vcf.gz -p NA12878 -r sim/CEU/simulated.art.hsxt.150l.100fc.400m.100s.NA12878.chr19.bam -m 60 -d reference/CEU.ped -o results/smartPhase.sim.NA12878.chr19.trio.AGV6UTR.allGeneRegionsCanonical.results.tsv -v -x -t