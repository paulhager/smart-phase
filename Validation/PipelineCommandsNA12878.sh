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
  $art -ss HSXt -i $f -p -l $length -f $cv -m $frag -s $sdev -o sim/YRI/simulated.art.hsxt.${length}l.${cv}fc.${frag}m.${sdev}s.${sample}. > ${sample}.out 2> ${sample}.err &
done

# Rename headers to prevent duplicate names
cat sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.1.1.fq | awk '{if(NR%4==1) $0="hapl1_"$0; print;}' > sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.1.1.fixedHeader.fq
cat sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.1.2.fq | awk '{if(NR%4==1) $0="hapl1_"$0; print;}' > sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.1.2.fixedHeader.fq

# Merge Haplotypes
cat sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr1.1.1.fq sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr1.2.1.fq > sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr1.merge.1.fq
cat sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr1.1.2.fq sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr1.2.2.fq > sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr1.merge.2.fq

cat sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.1.1.fixedHeader.fq sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.2.1.fq > sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.merge.1.fq
cat sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.1.2.fixedHeader.fq sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.2.2.fq > sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.merge.2.fq

# Map using BWA
bwa index reference/human_g1k_v37.fasta

bwa mem reference/human_g1k_v37.fasta sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr1.merge.1.fq sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr1.merge.2.fq > sim/CEU/simulated.art.NA12878.chr1.merge.pe.sam

bwa mem reference/human_g1k_v37.fasta sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.merge.1.fq sim/CEU/simulated.art.hsxt.150l.60fc.200m.10s.NA12878.chr19.merge.2.fq > sim/CEU/simulated.art.NA12878.chr19.merge.pe.sam

# SAM to BAM
samtools sort -o sim/CEU/simulated.art.NA12878.chr1.merge.pe.sort.bam sim/CEU/simulated.art.NA12878.chr1.merge.pe.sam

samtools sort -o sim/CEU/simulated.art.NA12878.chr19.merge.pe.sort.bam sim/CEU/simulated.art.NA12878.chr19.merge.pe.sam

# Extract those reads falling into u674
cut -c4- reference/u674control.rmd_DoC_minCov_8.chr1.bed > reference/u674control.rmd_DoC_minCov_8.chr1.cut.bed
samtools view -b -L reference/u674control.rmd_DoC_minCov_8.chr1.cut.bed sim/CEU/simulated.art.NA12878.chr1.merge.pe.sort.bam > sim/CEU/simulated.art.NA12878.chr1.u674.merge.pe.bam
samtools index sim/CEU/simulated.art.NA12878.chr1.u674.merge.pe.bam

cut -c4- reference/u674control.rmd_DoC_minCov_8.chr19.bed > reference/u674control.rmd_DoC_minCov_8.chr19.cut.bed
samtools view -b -L reference/u674control.rmd_DoC_minCov_8.chr19.cut.bed sim/CEU/simulated.art.NA12878.chr19.merge.pe.sort.bam > sim/CEU/simulated.art.NA12878.chr19.u674.merge.pe.bam
samtools index sim/CEU/simulated.art.NA12878.chr19.u674.merge.pe.bam

# Extract variants falling into GENCODE coding exons
bedtools intersect -a reference/u674control.rmd_DoC_minCov_8.chr1.bed -b reference/hg19_GENCODEv27lift37_codingExons_sorted_nochrM_merged.bed > reference/intersect.u674.hg19GENCODE_codingExons.chr1.bed
cut -c4- reference/intersect.u674.hg19GENCODE_codingExons.chr1.bed > reference/intersect.u674.hg19GENCODE_codingExons.chr1.cut.bed
vcftools --vcf sim/CEU/sim.CEU.trio.chr1.phased.vcf --bed reference/intersect.u674.hg19GENCODE_codingExons.chr1.cut.bed --out sim/CEU/sim.CEU.trio.chr1.phased.u674.hg19GENCODE_codingExons --recode --keep-INFO-all
bgzip -c sim/CEU/sim.CEU.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf > sim/CEU/sim.CEU.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz
tabix sim/CEU/sim.CEU.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz

bedtools intersect -a reference/u674control.rmd_DoC_minCov_8.chr19.bed -b reference/hg19_GENCODEv27lift37_codingExons_sorted_nochrM_merged.bed > reference/intersect.u674.hg19GENCODE_codingExons.chr19.bed
cut -c4- reference/intersect.u674.hg19GENCODE_codingExons.chr19.bed > reference/intersect.u674.hg19GENCODE_codingExons.chr19.cut.bed
vcftools --vcf sim/CEU/sim.CEU.trio.chr19.phased.vcf --bed reference/intersect.u674.hg19GENCODE_codingExons.chr19.cut.bed --out sim/CEU/sim.CEU.trio.chr19.phased.u674.hg19GENCODE_codingExons --recode --keep-INFO-all
bgzip -c sim/CEU/sim.CEU.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf > sim/CEU/sim.CEU.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz
tabix sim/CEU/sim.CEU.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz

# SmartPhase
java -jar smartPhase.jar -g=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/allGeneRegionsCanonical.HG19.GRCh37.chr1.bed -f=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz -a=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz -p=NA19240 -r=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/simulated.art.NA19240.chr1.u674.merge.pe.bam -m=60 -d=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/YRI.ped -o=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/smartPhase.sim.NA19240.chr1.trio.u674.GENCODE_cE.canonicalGenes.results.txt -x -h -v -t

java -jar smartPhase.jar -g=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/allGeneRegionsCanonical.HG19.GRCh37.chr19.bed -f=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz -a=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz -p=NA19240 -r=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/simulated.art.NA19240.chr19.u674.merge.pe.bam -m=60 -d=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/YRI.ped -o=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/smartPhase.sim.NA19240.chr19.trio.u674.GENCODE_cE.canonicalGenes.results.txt -x -h -v -t