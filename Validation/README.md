# Validation of SmartPhase on simulated data

The validation of SmartPhase on simulated data is based on two pipelines. 
The simulation pipeline generates FASTQ files by simulating reads of whole-exome sequencing separately for both haplotypes of the children of the CEU and the YRI trio.
The validation pipeline maps the generated reads, executes SmartPhase and WhatsHap and evaluates the phasing predictions.

Please note, that the simulation pipeline critically depends on the `V6_UTR_Probes_chr1_chr19.txt` file of the SureSelect Human All Exon V6 kit which must be requested from Agilent.
The resulting files of the simulation pipeline are available at [Google Drive](https://drive.google.com/drive/folders/1PLmow_1GPOi1enHDsOanQ3w2xVmUxlTh?usp=sharing) and can be directly used to run the validation pipeline.


## Requirements
- samtools, bedtools and vcftools and bwa (v0.7.15) must be on path
- Python 3.x including BioPython
- Java 10
- SHAPEIT:
	To download it go to <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download> and download version v2.r837 (shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz).
	Unpack the tgz file and place/link the `bin/shapeit binary` into `./bin/shapeit`.
- Wessim2: 
	Wessim will be cloned into the `Validation` directory from its [GitHub repository](https://github.com/sak042/Wessim) during step 1 of the pipeline.
	Make sure that all listed requirements are met. Further information can be found at <http://sak042.github.io/Wessim/>.
	Please note, in order to run Wessim2 as we did it in our pipeline you need the `_probes.txt` file of the underlying exome enrichment kit.
	We are using the SureSelect Human All Exon V6+UTR kit. Due to proprietary reasons we are not allowed to share the used `V6_UTR_Probes_chr1_chr19.txt` file. 
	Please contact [Agilent](https://www.agilent.com/en-us/contact-us/page) directly to ask for the file if you want to rerun the simulation pipeline completely.
	Otherwise you can simply use our generated FASTQ files provided at [Google Drive](https://drive.google.com/drive/folders/1PLmow_1GPOi1enHDsOanQ3w2xVmUxlTh?usp=sharing).
- GemSim:
	The error model files are required for NGS read generation. Please download GemSim from <http://sourceforge.net/projects/gemsim/> and save it in the `Validation` folder.
	The model file `ill100v5_p.gzip` is used in the pipeline.
- WhatsHap (only required if you want to run the comparison to WhatsHap):
	To install it visit <https://whatshap.readthedocs.io/en/latest/installation.html>.

## Simulation pipeline steps
1. Download reference data, Wessim repository and WhatsHap phasing comparison experiment scripts provided in the [phasing-comparison-experiments BitBucket repository of WhatsHap](https://bitbucket.org/whatshap/phasing-comparison-experiments/ "WhatsHap phasing comparison experiments").
2. Prepare a BED file representing WES captured protein-coding exons as the intersection of the regions captured by the Agilent SureSelect Human All Exon V6 kit and the regions of canonical protein-coding genes downloaded from the USCS table browser.
Both files are located in the [BED directory](https://github.com/paulhager/smart-phase/tree/master/BED).
3. Download and prepare the VCFs of the CEU and YRI trio provided by the [1000 Genomes FTP server](https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20140625_high_coverage_trios_broad/).
4. Run SHAPEIT for phasing both trios on chr1 and chr19.
5. Use WhatsHap script `artificial-child.py` to create an artificial child VCF based on the phased parents which is then used to generate an artifial trio VCF.
6. Use WhatsHap script `genomesimulator.py` to create haplotype FASTA files for the artifical child.
7. Extract variants of the artifical trio which are located in the regions of the BED file prepared in step 2.
8. Use the Wessim2 read simulation tool to create paired-end sequencing data with 100 bp read length and an average depth of coverage of 126.
9. Merge the FASTQ files generated for each haplotype.

## Validation pipeline steps
1. Run BWA-MEM to align simulated reads to the reference genome hg19.
2. Run SmartPhase for the children of the CEU and YRI trio in read-only and read-and-trio mode each on chr1 and chr19 using the BAM files produced in step 10 and the VCF files of step 11 and extract the metrics described in the manuscript.
3. Run WhatsHap for the children of the CEU and YRI trio in read-only and read-and-trio mode each on chr1 and chr19 using the BAM files produced in step 10 and the VCF files of step 11 and determine phasing errors.

You can download the simulated data used in the publication [here](http://ibis.helmholtz-muenchen.de/smartphase/smartphase_simulation_data.tar.gz).
