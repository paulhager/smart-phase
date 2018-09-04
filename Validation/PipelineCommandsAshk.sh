Pipeline Commands Ashk

@@@WHATSHAP@@@

DOWNLOAD LOCATIONS:
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_CallsIn2Technologies_05182015/HG002-multiall-fullcombine.vcf.gz
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_CallsIn2Technologies_05182015/HG003-multiall-fullcombine.vcf.gz
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_CallsIn2Technologies_05182015/HG004-multiall-fullcombine.vcf.gz
http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz
http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz
http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/genetic_map_chr1_combined_b37.txt
http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3.sample

STARTING CONFIGURATION

.//AshkenazimTrio
	.//AshkenazimTrio/father
		.//AshkenazimTrio/father/HG003-multiall-fullcombine.unphased.vcf.gz
		.//AshkenazimTrio/father/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
	.//AshkenazimTrio/mother
		.//AshkenazimTrio/mother/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz
		.//AshkenazimTrio/mother/HG004-multiall-fullcombine.unphased.vcf.gz
	.//AshkenazimTrio/son
		.//AshkenazimTrio/son/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
		.//AshkenazimTrio/son/HG002-multiall-fullcombine.unphased.vcf.gz
.//shapeit.v2.904.2.6.32-696.18.7.el6.x86_64
.//whatshap-comparison-experiments
.//reference
	..//1000GP_Phase3.sample
	..//genetic_map_chr1_combined_b37.txt
	..//1000GP_Phase3_chr1.legend.gz
	..//1000GP_Phase3_chr1.hap.gz
	..//AGV6UTR_covered_merged.bed


PIPELINE - TO BE EXECUTED FROM START DIRECTORY

# Unzip unphased vcfs
gunzip -c AshkenazimTrio/son/HG002-multiall-fullcombine.unphased.vcf.gz > AshkenazimTrio/son/HG002-multiall-fullcombine.unphased.vcf
gunzip -c AshkenazimTrio/father/HG003-multiall-fullcombine.unphased.vcf.gz > AshkenazimTrio/father/HG003-multiall-fullcombine.unphased.vcf
gunzip -c AshkenazimTrio/mother/HG004-multiall-fullcombine.unphased.vcf.gz > AshkenazimTrio/mother/HG004-multiall-fullcombine.unphased.vcf

# Merge vcfs using WhatsHap vcf_merge script
whatshap-comparison-experiments/scripts/vcf_merge_trio.pl AshkenazimTrio/mother/HG004-multiall-fullcombine.unphased.vcf AshkenazimTrio/father/HG003-multiall-fullcombine.unphased.vcf AshkenazimTrio/son/HG002-multiall-fullcombine.unphased.vcf > AshkenazimTrio/ashk.trio.unphased.vcf

# Extract chr1 from vcf
""awk '/^#/ || ($1 == "1")' AshkenazimTrio/ashk.trio.unphased.vcf > AshkenazimTrio/ashk.trio.unphased.chr1.vcf""

""awk '/^#/ || ($1 == "19")' AshkenazimTrio/ashk.trio.unphased.vcf > AshkenazimTrio/ashk.trio.unphased.chr19.vcf""

# Changes spaces to tabs because the merge script kills all tabs
awk 'NR>315 {gsub(" ","\t")}1' AshkenazimTrio/ashk.trio.unphased.chr1.vcf > AshkenazimTrio/ashk.trio.unphased.chr1.tab.vcf

awk 'NR>315 {gsub(" ","\t")}1' AshkenazimTrio/ashk.trio.unphased.chr19.vcf > AshkenazimTrio/ashk.trio.unphased.chr19.tab.vcf

# Tabix sort chr1 vcf file
bgzip -c AshkenazimTrio/ashk.trio.unphased.chr1.tab.vcf > AshkenazimTrio/ashk.trio.unphased.chr1.tab.vcf.gz
tabix AshkenazimTrio/ashk.trio.unphased.chr1.tab.vcf.gz

bgzip -c AshkenazimTrio/ashk.trio.unphased.chr19.tab.vcf > AshkenazimTrio/ashk.trio.unphased.chr19.tab.vcf.gz
tabix AshkenazimTrio/ashk.trio.unphased.chr19.tab.vcf.gz

# Run SHAPEIT Check, then SHAPEIT, then convert to vcf
shapeit -check -V AshkenazimTrio/ashk.trio.unphased.chr1.tab.vcf -M reference/genetic_map_chr1_combined_b37.txt --input-ref reference/1000GP_Phase3_chr1.hap.gz reference/1000GP_Phase3_chr1.legend.gz reference/1000GP_Phase3.sample --output-log shapeit/trio.chr1 || true) > shapeit/trio.chr1.check.log 2>&1
shapeit -V AshkenazimTrio/ashk.trio.unphased.chr1.tab.vcf --exclude-snp shapeit/trio.chr1.snp.strand.exclude' -M reference/genetic_map_chr1_combined_b37.txt --input-ref reference/1000GP_Phase3_chr1.hap.gz reference/1000GP_Phase3_chr1.legend.gz reference/1000GP_Phase3.sample -O shapeit/trio.chr1 > shapeit/trio.chr1.run.log 2>&
shapeit -convert --input-haps shapeit/trio.chr1 --output-vcf shapeit/trio.chr1.phased.vcf > shapeit/trio.chr1.phased.vcf.log 2>&1

shapeit -check -V AshkenazimTrio/ashk.trio.unphased.chr19.tab.vcf -M reference/genetic_map_chr1_combined_b37.txt --input-ref reference/1000GP_Phase3_chr19.hap.gz reference/1000GP_Phase3_chr19.legend.gz reference/1000GP_Phase3.sample --output-log shapeit/trio.chr19 || true) > shapeit/trio.chr19.check.log 2>&1
shapeit -V AshkenazimTrio/ashk.trio.unphased.chr19.tab.vcf --exclude-snp shapeit/trio.chchr19r1.snp.strand.exclude' -M reference/genetic_map_chr19_combined_b37.txt --input-ref reference/1000GP_Phase3_chr19.hap.gz reference/1000GP_Phase3_chr19.legend.gz reference/1000GP_Phase3.sample -O shapeit/trio.chr19 > shapeit/trio.chr19.run.log 2>&
shapeit -convert --input-haps shapeit/trio.chr19 --output-vcf shapeit/trio.chr19.phased.vcf > shapeit/trio.chr19.phased.vcf.log 2>&1

# Split VCF into three
vcf-subset -c HG002 ashk.trio.unphased.chr1.tab_phased.vcf > ashk.individual.child.chr1.tab.phased.vcf
vcf-subset -c HG003 ashk.trio.unphased.chr1.tab_phased.vcf > ashk.individual.father.chr1.tab.phased.vcf
vcf-subset -c HG004 ashk.trio.unphased.chr1.tab_phased.vcf > ashk.individual.mother.chr1.tab.phased.vcf

bcftools view -s HG002 AshkenazimTrio/ashk.trio.unphased.chr19.tab_biallelic_phased.vcf.gz > AshkenazimTrio/ashk.individual.child.chr19.tab.phased.vcf
bcftools view -s HG003 AshkenazimTrio/ashk.trio.unphased.chr19.tab_biallelic_phased.vcf.gz > AshkenazimTrio/ashk.individual.father.chr19.tab.phased.vcf
bcftools view -s HG004 AshkenazimTrio/ashk.trio.unphased.chr19.tab_biallelic_phased.vcf.gz > AshkenazimTrio/ashk.individual.mother.chr19.tab.phased.vcf

# Multiply distances by factor (What is the point of this?)
awk 'NR==1 {{print}} NR>1 {{print $1, $2*1.5, $3*1.5 }}' reference/genetic_map_chr1_combined_b37.txt > reference/genetic_map_chr1.txt

# Create symbolic links
ln -s shapeit/ashk.individual.mother.chr1.tab.phased.vcf sim.individual.mother.chr1.tab.phased.vcf
ln -s shapeit/ashk.individual.father.chr1.tab.phased.vcf sim.individual.father.chr1.tab.phased.vcf

ln -s shapeit/ashk.individual.mother.chr19.tab.phased.vcf AshkenazimTrio/sim.individual.mother.chr19.tab.phased.vcf
ln -s shapeit/ashk.individual.father.chr19.tab.phased.vcf AshkenazimTrio/sim.individual.father.chr19.tab.phased.vcf

# Create artificial child
whatshap-comparison-experiments/scripts/artificial-child.py reference/genetic_map_chr1.txt shapeit/sim.individual.mother.chr1.tab.phased.vcf shapeit/sim.individual.father.chr1.tab.phased.vcf HG002 sim/mother.chr1.true.recomb sim/father.chr1.true.recomb > vcf/sim.child.chr1.phased.vcf 2>log

whatshap-comparison-experiments/scripts/artificial-child.py reference/genetic_map_chr19.txt shapeit/sim.individual.mother.chr19.tab.phased.vcf shapeit/sim.individual.father.chr19.tab.phased.vcf HG002 sim/mother.chr19.true.recomb sim/father.chr19.true.recomb > vcf/sim.child.chr19.phased.vcf 2>log

# Merge artifical trio
bgzip -c shapeit/sim.individual.mother.chr1.tab.phased.vcf > vcf/sim.individual.mother.chr1.tab.phased.vcf.gz
bgzip -c shapeit/sim.individual.father.chr1.tab.phased.vcf > vcf/sim.individual.father.chr1.tab.phased.vcf.gz
bgzip -c vcf/sim.child.chr1.phased.vcf > vcf/sim.child.chr1.phased.vcf.gz
tabix vcf/sim.individual.mother.chr1.tab.phased.vcf.gz 
tabix vcf/sim.individual.father.chr1.tab.phased.vcf.gz
tabix vcf/sim.child.chr1.phased.vcf.gz 
vcf-merge vcf/sim.individual.mother.chr1.tab.phased.vcf.gz vcf/sim.individual.father.chr1.tab.phased.vcf.gz vcf/sim.child.chr1.phased.vcf.gz > vcf/sim.trio.chr1.phased.vcf

bgzip -c AshkenazimTrio/ashk.individual.mother.chr19.tab.phased.vcf > vcf/sim.individual.mother.chr19.tab.phased.vcf.gz
bgzip -c AshkenazimTrio/ashk.individual.father.chr19.tab.phased.vcf > vcf/sim.individual.father.chr19.tab.phased.vcf.gz
bgzip -c AshkenazimTrio/ashk.individual.child.chr19.tab.phased.vcf > vcf/sim.individual.child.chr19.tab.phased.vcf.gz
tabix vcf/sim.individual.mother.chr19.tab.phased.vcf.gz 
tabix vcf/sim.individual.father.chr19.tab.phased.vcf.gz
tabix vcf/sim.individual.child.chr19.tab.phased.vcf.gz 
vcf-merge vcf/sim.individual.mother.chr19.tab.phased.vcf.gz vcf/sim.individual.father.chr19.tab.phased.vcf.gz vcf/sim.individual.child.chr19.tab.phased.vcf.gz > vcf/sim.trio.chr19.phased.vcf

# Create child true haplotype fastas
mkdir -p sim/tmp

whatshap-comparison-experiments/scripts/genomesimulator.py -c 1 vcf/sim.child.chr1.phased.vcf reference/human_g1k_v37.fasta sim/tmp

whatshap-comparison-experiments/scripts/genomesimulator.py -c 19 AshkenazimTrio/ashk.individual.child.chr19.tab.phased.vcf reference/human_g1k_v37.fasta sim/tmp

# Create ped
echo family	HG002	HG003	HG004	0	2 > reference/ashk.ped

# Cut desired regions out of fasta
for f in sim/tmp/*fasta; do
  outFile=${f/.fasta/_kit.fasta}
  outFile=${outFile/tmp/HG002}
  bedtools getfasta -fi $f -bed reference/AGV6UTR_covered_merged.bed -fo $outFile
  rm $f
	rm $f.fai
done

# Simulate reads of both haplotypes using ART
art_bin_MountRainier/art_illumina -ss HSXt -i sim/HG002/HG002.chr1.1.fasta -p -l 150 -f 60 -m 200 -s 10 -o sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.1.
art_bin_MountRainier/art_illumina -ss HSXt -i sim/HG002/HG002.chr1.2.fasta -p -l 150 -f 60 -m 200 -s 10 -o sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.2.

art_bin_MountRainier/art_illumina -ss HSXt -i sim/HG002/HG002.chr19.1.fasta -p -l 150 -f 60 -m 200 -s 10 -o sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.1.
art_bin_MountRainier/art_illumina -ss HSXt -i sim/HG002/HG002.chr19.2.fasta -p -l 150 -f 60 -m 200 -s 10 -o sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.2.

# Rename headers to prevent duplicate names
cat sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.1.1.fq | awk '{if(NR%4==1) $0=$0"_hapl1"; print;}' > sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.1.1.fixedHeader.fq
cat sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.1.2.fq | awk '{if(NR%4==1) $0=$0"_hapl1"; print;}' > sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.1.2.fixedHeader.fq

# Merge Haplotypes
cat sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.1.1.fq sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.2.1.fq > sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.merge.1.fq
cat sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.1.2.fq sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.2.2.fq > sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.merge.2.fq

cat sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.1.1.fixedHeader.fq sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.2.1.fq > sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.merge.1.fq
cat sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.1.2.fixedHeader.fq sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.2.2.fq > sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.merge.2.fq

# Map using BWA
bwa index reference/human_g1k_v37.chr1.fasta

bwa mem reference/human_g1k_v37.fasta sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.merge.1.fq sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr1.merge.2.fq > sim/HG002/simulated.art.HG002.chr1.merge.pe.sam

bwa mem reference/human_g1k_v37.fasta sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.merge.1.fq sim/HG002/simulated.art.hsxt.150l.60fc.200m.10s.HG002.chr19.merge.2.fq > sim/HG002/simulated.art.HG002.chr19.merge.pe.sam

# SAM to BAM
samtools sort -o sim/HG002/simulated.art.HG002.chr1.merge.pe.sorted.bam sim/HG002/simulated.art.HG002.chr1.merge.pe.sam 

samtools sort -o sim/HG002/simulated.art.HG002.chr19.merge.pe.sorted.bam sim/HG002/simulated.art.HG002.chr19.merge.pe.sam 

# Extract reads in u674 region
cut -c4- reference/u674control.rmd_DoC_minCov_8.chr1.bed > reference/u674control.rmd_DoC_minCov_8.chr1.cut.bed
samtools view -b -L reference/u674control.rmd_DoC_minCov_8.chr1.cut.bed sim/HG002/simulated.art.HG002.chr1.merge.pe.sorted.bam > sim/HG002/simulated.art.HG002.chr1.u674control.merge.pe.sorted.bam
samtools index sim/HG002/simulated.art.HG002.chr1.u674control.merge.pe.sorted.bam

cut -c4- reference/u674control.rmd_DoC_minCov_8.chr19.bed > reference/u674control.rmd_DoC_minCov_8.chr19.cut.bed
samtools view -b -L reference/u674control.rmd_DoC_minCov_8.chr19.cut.bed sim/HG002/simulated.art.HG002.chr19.merge.pe.sorted.bam > sim/HG002/simulated.art.HG002.chr19.u674control.merge.pe.sorted.bam
samtools index sim/HG002/simulated.art.HG002.chr19.u674control.merge.pe.sorted.bam

# Extract variants in GENCODE coding exons region
bedtools intersect -a reference/u674control.rmd_DoC_minCov_8.chr1.bed -b reference/hg19_GENCODEv27lift37_codingExons_sorted_nochrM_merged.bed > reference/intersect.u674.hg19GENCODE_codingExons.chr1.bed
cut -c4- reference/intersect.u674.hg19GENCODE_codingExons.chr1.bed > reference/intersect.u674.hg19GENCODE_codingExons.chr1.cut.bed
vcftools --vcf vcf/sim.trio.chr1.phased.vcf --bed reference/intersect.u674.hg19GENCODE_codingExons.chr1.cut.bed --out vcf/sim.ASHK.trio.chr1.phased.u674.hg19GENCODE_codingExons --recode --keep-INFO-all
bgzip -c vcf/sim.ASHK.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf > vcf/sim.ASHK.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz
tabix vcf/sim.ASHK.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz

bedtools intersect -a reference/u674control.rmd_DoC_minCov_8.chr19.bed -b reference/hg19_GENCODEv27lift37_codingExons_sorted_nochrM_merged.bed > reference/intersect.u674.hg19GENCODE_codingExons.chr19.bed
cut -c4- reference/intersect.u674.hg19GENCODE_codingExons.chr19.bed > reference/intersect.u674.hg19GENCODE_codingExons.chr19.cut.bed
vcftools --vcf vcf/sim.trio.chr19.phased.vcf --bed reference/intersect.u674.hg19GENCODE_codingExons.chr19.cut.bed --out vcf/sim.ASHK.trio.chr19.phased.u674.hg19GENCODE_codingExons --recode --keep-INFO-all
bgzip -c vcf/sim.ASHK.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf > vcf/sim.ASHK.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz
tabix vcf/sim.ASHK.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz

# Run smartphase

java -jar smartPhase.jar -g=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/allGeneRegionsCanonical.HG19.GRCh37.chr1.bed -f=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.ASHK.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz -a=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.ASHK.trio.chr1.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz -p=HG002 -r=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/HG002/simulated.art.HG002.chr1.u674control.merge.pe.sorted.bam -m=60 -d=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/ashk.ped -o=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/smartPhase.sim.HG002.chr1.trio.u674.GENCODE_cE.canonicalGenes.results.txt -x -h -v -t

java -jar smartPhase.jar -g=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/allGeneRegionsCanonical.HG19.GRCh37.chr19.bed -f=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.ASHK.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz -a=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.ASHK.trio.chr19.phased.u674.hg19GENCODE_codingExons.recode.vcf.gz -p=HG002 -r=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/HG002/simulated.art.HG002.chr19.u674control.merge.pe.sorted.bam -m=60 -d=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/ashk.ped -o=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/smartPhase.sim.HG002.chr19.trio.u674.GENCODE_cE.canonicalGenes.results.txt -x -h -v -t

#OLD

cut -c4- reference/AGV6UTR_covered_merged.bed > reference/AGV6UTR_covered_merged.chrCUT.bed
samtools view -b -L reference/AGV6UTR_covered_merged.chrCUT.bed sim/simulated.art.HG002.chr1.merge.pe.sorted.bam > sim/simulated.art.HG002.chr1.AGV6UTR.merge.pe.sorted.bam

cut -c4- reference/allGeneRegionsCanonical.HG19.GRCh37.bed > reference/allGeneRegionsCanonical.HG19.GRCh37.chrCut.bed
samtools view -b -L reference/allGeneRegionsCanonical.HG19.GRCh37.chrCUT.bed sim/simulated.art.HG002.chr1.merge.pe.sorted.bam > sim/simulated.art.HG002.chr1.CanonicalGeneRegions.merge.pe.sorted.bam

cut -c4- reference/UCSCCanonicalExons.chr1.bed > reference/UCSCCanonicalExons.chr1.cut.bed
samtools view -b -L reference/UCSCCanonicalExons.chr1.cut.bed sim/simulated.art.HG002.chr1.merge.pe.sorted.bam > sim/simulated.art.HG002.chr1.CanonicalExons.merge.pe.sorted.bam


samtools index sim/simulated.art.HG002.chr1.AGV6UTR.merge.pe.sorted.bam 
samtools index sim/simulated.art.HG002.chr1.CanonicalGeneRegions.merge.pe.sorted.bam
samtools index sim/simulated.art.HG002.chr1.CanonicalExons.merge.pe.sorted.bam

bgzip -c vcf/sim.trio.chr1.phased.vcf > vcf/sim.trio.chr1.phased.vcf.gz
tabix vcf/sim.trio.chr1.phased.vcf.gz 

java -jar smartPhase.jar -g=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/AGV6UTR_covered_merged.chr1.bed -f=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.trio.chr1.phased.vcf.gz -a=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.trio.chr1.phased.vcf.gz -p=HG002 -r=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/simulated.art.HG002.chr1.AGV6UTR.merge.pe.sorted.bam -m=60  -t -d=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/ashk.ped -o=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/smartPhase.sim.HG002.chr1.trio.AGV6UTR.results.txt -x -h

java -jar smartPhase.jar -g=/Users/paulhager/Documents/SmartPhase/allGeneRegionsCanonical.HG19.GRCh37.chr1.bed -f=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.trio.chr1.phased.vcf.gz -a=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.trio.chr1.phased.vcf.gz -p=HG002 -r=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/simulated.art.HG002.chr1.CanonicalGeneRegions.merge.pe.sorted.bam -m=60  -t -d=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/ashk.ped -o=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/smartPhase.sim.HG002.chr1.trio.AllCanonicalGeneRegionsHG19GRCH37.results.txt -x -h

java -jar smartPhase.jar -g=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/UCSCCanonicalExons.chr1.bed -f=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.trio.chr1.phased.vcf.gz -a=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.trio.chr1.phased.vcf.gz -p=HG002 -r=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/simulated.art.HG002.chr1.CanonicalExons.merge.pe.sorted.bam -m=60  -t -d=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/ashk.ped -o=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/smartPhase.sim.HG002.chr1.trio.CanonicalExons.results.txt -x -h

java -jar smartPhase.jar -g=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/hg19_GENCODEv27lift37_codingExons_sorted_nochrM_merged.bed -f=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.trio.chr1.phased.vcf.gz -a=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/vcf/sim.trio.chr1.phased.vcf.gz -p=HG002 -r=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/sim/simulated.art.HG002.chr1.u674control.merge.pe.sorted.bam -m=60  -t -d=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/reference/ashk.ped -o=/Users/paulhager/Documents/SmartPhase/Publication/Comparison/results/smartPhase.sim.HG002.chr1.trio.u674.results.txt -x -h

