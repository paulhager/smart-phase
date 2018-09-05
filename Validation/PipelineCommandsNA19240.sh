Pipeline Commands NA19240

DOWNLOAD LOCATIONS:
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/

STARTING CONFIGURATION

NA19240
    YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.vcf.gz
    

# Unzip unphased vcfs
gunzip -c YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.vcf.gz > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.vcf

# Extract chr1 and chr19 from vcf
""awk '/^#/ || ($1 == "1")' YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.vcf > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1.vcf""

""awk '/^#/ || ($1 == "19")' YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.vcf > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19.vcf""

# Zip and index
bgzip -c YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1.vcf > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1.vcf.gz
tabix YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1.vcf.gz

bgzip -c YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19.vcf > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19.vcf.gz
tabix YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19.vcf.gz



# Run SHAPEIT Check, then SHAPEIT, then convert to vcf
shapeit -check -V AshkenazimTrio/ashk.trio.unphased.chr1.tab.vcf -M reference/genetic_map_chr1_combined_b37.txt --input-ref reference/1000GP_Phase3_chr1.hap.gz reference/1000GP_Phase3_chr1.legend.gz reference/1000GP_Phase3.sample --output-log shapeit/trio.chr1 || true) > shapeit/trio.chr1.check.log 2>&1
shapeit -V AshkenazimTrio/ashk.trio.unphased.chr1.tab.vcf --exclude-snp shapeit/trio.chr1.snp.strand.exclude' -M reference/genetic_map_chr1_combined_b37.txt --input-ref reference/1000GP_Phase3_chr1.hap.gz reference/1000GP_Phase3_chr1.legend.gz reference/1000GP_Phase3.sample -O shapeit/trio.chr1 > shapeit/trio.chr1.run.log 2>&
shapeit -convert --input-haps shapeit/trio.chr1 --output-vcf YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1_biallelic_phased.vcf > shapeit/trio.chr1.phased.vcf.log 2>&1

shapeit -check -V AshkenazimTrio/ashk.trio.unphased.chr19.tab.vcf -M reference/genetic_map_chr19_combined_b37.txt --input-ref reference/1000GP_Phase3_chr19.hap.gz reference/1000GP_Phase3_chr19.legend.gz reference/1000GP_Phase3.sample --output-log shapeit/trio.chr19 || true) > shapeit/trio.chr19.check.log 2>&1
shapeit -V AshkenazimTrio/ashk.trio.unphased.chr19.tab.vcf --exclude-snp shapeit/trio.chr19.snp.strand.exclude' -M reference/genetic_map_chr19_combined_b37.txt --input-ref reference/1000GP_Phase3_chr19.hap.gz reference/1000GP_Phase3_chr19.legend.gz reference/1000GP_Phase3.sample -O shapeit/trio.chr19 > shapeit/trio.chr19.run.log 2>&
shapeit -convert --input-haps shapeit/trio.chr19 --output-vcf YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19_biallelic_phased.vcf > shapeit/trio.chr19.phased.vcf.log 2>&1

# Split VCF into three
bcftools view -s NA19240 YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1_biallelic_phased.vcf.gz > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19240.chr1_biallelic_phased.vcf
bcftools view -s NA19238 YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1_biallelic_phased.vcf.gz > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr1_biallelic_phased.vcf
bcftools view -s NA19239 YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr1_biallelic_phased.vcf.gz > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr1_biallelic_phased.vcf

bcftools view -s NA19240 YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19_biallelic_phased.vcf.gz > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19240.chr19_biallelic_phased.vcf
bcftools view -s NA19238 YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19_biallelic_phased.vcf.gz > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr19_biallelic_phased.vcf
bcftools view -s NA19239 YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.chr19_biallelic_phased.vcf.gz > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr19_biallelic_phased.vcf

# Multiply distances by factor (What is the point of this?)
awk 'NR==1 {{print}} NR>1 {{print $1, $2*1.5, $3*1.5 }}' reference/genetic_map_chr1_combined_b37.txt > reference/genetic_map_chr1.txt

awk 'NR==1 {{print}} NR>1 {{print $1, $2*1.5, $3*1.5 }}' reference/genetic_map_chr19_combined_b37.txt > reference/genetic_map_chr19.txt

# Create artifical child
whatshap-comparison-experiments/scripts/artificial-child.py reference/genetic_map_chr1.txt YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr1_biallelic_phased.vcf YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr1_biallelic_phased.vcf NA19240 sim/YRI/YRI.NA19239.chr1.true.recomb sim/YRI/YRI.NA19238.chr1.true.recomb > sim/YRI/sim.YRI.NA19240.chr1.phased.vcf

whatshap-comparison-experiments/scripts/artificial-child.py reference/genetic_map_chr19.txt YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr19_biallelic_phased.vcf YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr19_biallelic_phased.vcf NA19240 sim/YRI/YRI.NA19239.chr19.true.recomb sim/YRI/YRI.NA19238.chr19.true.recomb > sim/YRI/sim.YRI.NA19240.chr19.phased.vcf

# Merge artifical trio
bgzip -c YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr1_biallelic_phased.vcf > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr1_biallelic_phased.vcf.gz
tabix YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr1_biallelic_phased.vcf.gz
bgzip -c YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr1_biallelic_phased.vcf > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr1_biallelic_phased.vcf.gz
tabix YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr1_biallelic_phased.vcf.gz
bgzip -c sim/YRI/sim.YRI.NA19240.chr1.phased.vcf > sim/YRI/sim.YRI.NA19240.chr1.phased.vcf.gz
tabix sim/YRI/sim.YRI.NA19240.chr1.phased.vcf.gz

vcf-merge YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr1_biallelic_phased.vcf.gz YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr1_biallelic_phased.vcf.gz sim/YRI/sim.YRI.NA19240.chr1.phased.vcf.gz > sim/YRI/sim.YRI.trio.chr1.phased.vcf


bgzip -c YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr19_biallelic_phased.vcf > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr19_biallelic_phased.vcf.gz
tabix YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr19_biallelic_phased.vcf.gz
bgzip -c YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr19_biallelic_phased.vcf > YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr19_biallelic_phased.vcf.gz
tabix YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr19_biallelic_phased.vcf.gz
bgzip -c sim/YRI/sim.YRI.NA19240.chr19.phased.vcf > sim/YRI/sim.YRI.NA19240.chr19.phased.vcf.gz
tabix sim/YRI/sim.YRI.NA19240.chr19.phased.vcf.gz

vcf-merge YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19239.chr19_biallelic_phased.vcf.gz YRI/YRI.wgs.consensus.20131118.snps_indels.high_coverage_pcr_free.genotypes.NA19238.chr19_biallelic_phased.vcf.gz sim/YRI/sim.YRI.NA19240.chr19.phased.vcf.gz > sim/YRI/sim.YRI.trio.chr19.phased.vcf

# Create child true haplotype fastas
mkdir -p sim/tmp

whatshap-comparison-experiments/scripts/genomesimulator.py -c 1 sim/YRI/sim.YRI.NA19240.chr1.phased.vcf reference/human_g1k_v37.fasta sim/tmp/

whatshap-comparison-experiments/scripts/genomesimulator.py -c 19 sim/YRI/sim.YRI.NA19240.chr19.phased.vcf reference/human_g1k_v37.fasta sim/tmp/

# Create PED
echo YRI	NA19240	NA19238	NA19239	0	2 > reference/YRI.ped

# Cut desired regions out of fasta
for f in sim/tmp/*fasta; do
  outFile=${f/.fasta/_kit.fasta}
  outFile=${outFile/tmp/YRI}
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

for f in sim/YRI/*_kit.fasta; do
  sample=$( basename $f | cut -d "_" -f1 )
  echo start for $sample
  $art -ss HSXt -i $f -p -l $length -f $cv -m $frag -s $sdev -o sim/YRI/simulated.art.hsxt.${length}l.${cv}fc.${frag}m.${sdev}s.${sample}. > ${sample}.out 2> ${sample}.err &
done

# Rename headers to prevent duplicate names
for f in sim/YRI/*.1.[12].fq; do
  out=${f/.fq/.fixed_header.fq.gz}
  echo $f to $out
  cat $f | awk '{if(NR%4==1) $0="@hapl1_"$0; print;}' | gzip > $out &
done

for f in sim/YRI/*.2.[12].fq; do
  out=${f/.fq/.fixed_header.fq.gz}
  echo $f to $out
  cat $f | awk '{if(NR%4==1) $0="@hapl2_"$0; print;}' | gzip > $out &
done

# Merge Haplotypes
for f in sim/YRI/*1.1.fixed_header.fq.gz; do
  f2=${f/1.1.fixed_header.fq.gz/2.1.fixed_header.fq.gz}
  gzcat $f $f2 | gzip > ${f/1.1.fixed_header.fq.gz/1.fq.gz} &
done

for f in sim/YRI/*1.2.fixed_header.fq.gz; do
  f2=${f/1.2.fixed_header.fq.gz/2.2.fixed_header.fq.gz}
  gzcat $f $f2 | gzip > ${f/1.2.fixed_header.fq.gz/2.fq.gz} &
done

# Map using BWA
bwa index reference/human_g1k_v37.chr1.fasta

for f in sim/YRI/*1.fq.gz; do
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

# Extract variants in kit and intersection of kit with gencode coding exons regions
bedtools intersect -a reference/AGV6UTR_covered_merged.bed -b reference/hg19_GENCODEv27lift37_codingExons_sorted_nochrM_merged.bed > reference/intersect.AGV6UTR.hg19GENCODE_codingExons.bed

cut -c4- reference/intersect.AGV6UTR.hg19GENCODE_codingExons.bed > reference/intersect.AGV6UTR.hg19GENCODE_codingExons.cut.bed
cut -c4- reference/AGV6UTR_covered_merged.bed > reference/AGV6UTR_covered_merged.cut.bed

for f in sim/YRI/sim.YRI.trio.*.phased.vcf; do
  out1=${f/.vcf/.AGV6UTR.hg19GENCODE_codingExons}
  out2=${f/.vcf/.AGV6UTR}

  vcftools --vcf $f --bed reference/intersect.AGV6UTR.hg19GENCODE_codingExons.cut.bed --out $out1 --recode --keep-INFO-all
  vcftools --vcf $f --bed reference/AGV6UTR_covered_merged.cut.bed --out $out2 --recode --keep-INFO-all

  out1=$out1.recode.vcf
  out2=$out2.recode.vcf

  bgzip -c $out1 > $out1.gz
  bgzip -c $out2 > $out2.gz

  tabix $out1.gz
  tabix $out2.gz
done

# Run SmartPhase

java -jar smartPhase.jar -g /Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/allGeneRegionsCanonical.HG19.GRCh37.bed -f /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr1.phased.AGV6UTR.hg19GENCODE_codingExons.recode.vcf.gz -a /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr1.phased.AGV6UTR.recode.vcf.gz -p NA19240 -r /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/simulated.art.hsxt.150l.100fc.400m.100s.NA19240.chr1.bam -m 60 -d /Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/YRI.ped -o /Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/YRI/smartPhase.sim.NA19240.chr1.trio.AGV6UTR.GENCODE_cE.canonicalGenes.results.tsv -x -v -t

java -jar smartPhase.jar -g /Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/allGeneRegionsCanonical.HG19.GRCh37.bed -f /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr1.phased.AGV6UTR.hg19GENCODE_codingExons.recode.vcf.gz -a /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr1.phased.AGV6UTR.recode.vcf.gz -p NA19240 -r /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/simulated.art.hsxt.150l.100fc.400m.100s.NA19240.chr1.bam -m 60 -d /Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/YRI.ped -o /Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/YRI/smartPhase.sim.NA19240.chr1.NOtrio.AGV6UTR.GENCODE_cE.canonicalGenes.results.tsv -x -v

java -jar smartPhase.jar -g /Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/allGeneRegionsCanonical.HG19.GRCh37.bed -f /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr19.phased.AGV6UTR.hg19GENCODE_codingExons.recode.vcf.gz -a /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr19.phased.AGV6UTR.recode.vcf.gz -p NA19240 -r /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/simulated.art.hsxt.150l.100fc.400m.100s.NA19240.chr19.bam -m 60 -d /Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/YRI.ped -o /Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/YRI/smartPhase.sim.NA19240.chr19.trio.AGV6UTR.GENCODE_cE.canonicalGenes.results.tsv -x -v -t

java -jar smartPhase.jar -g /Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/allGeneRegionsCanonical.HG19.GRCh37.bed -f /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr19.phased.AGV6UTR.hg19GENCODE_codingExons.recode.vcf.gz -a /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/sim.YRI.trio.chr19.phased.AGV6UTR.recode.vcf.gz -p NA19240 -r /Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/YRI/simulated.art.hsxt.150l.100fc.400m.100s.NA19240.chr19.bam -m 60 -d /Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/YRI.ped -o /Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/YRI/smartPhase.sim.NA19240.chr19.NOtrio.AGV6UTR.GENCODE_cE.canonicalGenes.results.tsv -x -v
