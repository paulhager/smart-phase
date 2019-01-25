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
