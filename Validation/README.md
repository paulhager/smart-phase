# Validation of SmartPhase on simulated data

## Requirements
- samtools, bedtools and vcftools and bwa (v0.7.15) must be on path
- Python 3.x including BioPython
- SHAPEIT:
	To download it go to <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download> and download version v2.r837 (shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz).
	Unpack the tgz file and place/link the `bin/shapeit binary` into `./bin/shapeit`.
- ART simulation tools: 
	To download it go to <https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm> and download version ART-MountRainier-2016-06-05 (ART-bin-MountRainier-2016.06.05-Linux64.tgz).
	Unpack the tgz file and place/link the `art_bin_MountRainier/art_illumina` into `./bin/art_illumina`.
- Java 10
- WhatsHap (only required if you want to run the comparison to WhatsHap):
	To install it visit <https://whatshap.readthedocs.io/en/latest/installation.html>.

## Pipeline steps
1. Download reference data and scripts provided in the [phasing-comparison-experiments BitBucket repository of WhatsHap](https://bitbucket.org/whatshap/phasing-comparison-experiments/ "WhatsHap phasing comparison experiments").
2. Preparing a BED file representing WES captured protein-coding exons as the intersection of the regions captured by the Agilent SureSelect Human All Exon V6 kit and the regions of canonical protein-coding genes downloaded from the USCS table browser.
Both files are located in the [BED directory](https://github.com/paulhager/smart-phase/tree/master/BED).
3. Downloading the VCFs for the CEU and YRI trio provided by the [1000 Genomes FTP server](https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/).
4. Run SHAPEIT for phasing both trios on chr1 and chr19.
5. Use WhatsHap script `artificial-child.py` to create an artificial child VCF based on the phased parents which is then used to generate an artifial trio VCF.
6. Use WhatsHap script `genomesimulator.py` to create haplotype FASTA files for the artifical child.
7. Cut FASTA files to only include regions captured by the Agilent SureSelect Human All Exon V6 kit.
8. Use the ART Illumina read simulation tool to create paired-end sequencing data with 150 bp read length and an average depth of coverage of 100.
9. Merge the FASTQ files generated for each haplotype.
10. Run BWA-MEM to align simulated reads to the reference genome hg19.
11. Extract variants of the artifical trio which are located in the regions of the BED file prepared in step 2.
12. Run SmartPhase for the children of the CEU and YRI trio in read-only and read-and-trio mode each on chr1 and chr19 using the BAM files produced in step 10 and the VCF files of step 11 and extract the metrics described in the manuscript.
13. Run WhatsHap for the children of the CEU and YRI trio in read-only and read-and-trio mode each on chr1 and chr19 using the BAM files produced in step 10 and the VCF files of step 11 and determine phasing errors.

You can download the simulated data used in the publication [here](http://ibis.helmholtz-muenchen.de/smartphase/smartphase_simulation_data.tar.gz).
