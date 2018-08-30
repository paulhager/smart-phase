# smart-phase

SmartPhase is a phasing tool tailored for clinical use in genetic diagnosis pipelines. It accurately and efficiently reduces the number of possible compound heterozygous variant pairs being examined around either genetic loci or from a list of preselected variant pairs. To achieve this, SmartPhase is able to incorperate parental genotype information as well as reads generated from either DNA- or RNA-sequencing. Furthermore, it incorporates existing haplotype information and applies logical rules to exclude variant constellations that cannot be disease causing.

For a more thorough explanation of SmartPhase and it's validation, please refer to the following: <LINK TO ARTICLE/PREPRINT>

SmartPhase is offered as an executable jar in addition to its source code. The following serves as a general overview and explanation of all its arguments and recommended values.

# Requirements

- Java 8

# Usage

SmartPhase can either be run in explorative or analytic mode. In explorative mode, SmartPhase parses provided genomic regions of interest (typically protein coding regions) provided in a bed file and phases all variants of the filtered variants file that fall within a region to each other. In analytic mode, no regions of interest must be given, just specific filtered variants that should be phased to one another and the regions are created around them. Analytic mode lends itself more to an analysis of a cohort where heavy pre-filtering has already been done and only specific candidate variant pairs should be analysed whereas explorative mode is more all-purpose and can be run at any point in the analysis pipeline.

Basic usage:

java -jar SmartPhase.jar -g /path/to/genomic/regions/to/be/phased/file/allGeneRegionsCanonical.HG19.GRCh37.bed -a /path/to/vcf/containing/all/variants/sample.vcf.gz -p PID12345 -r /path/to/bams/containing/reads/DNAseq.reads.bam,/path/to/bams/containing/reads/RNAseq.reads.bam -m 60,255 -d /path/to/ped/family.ped -o /path/to/desired/output/file/output.tsv -t

Options:

-f or --filtered-variants

The full path to the file containing those variants that should be phased. If running in analytic mode, this file is required. If running in explorative mode, this can be left blank and will be set equal to the all-variants file provided, thus assuming that all variants that fall within the genomic regions of interest should be phased. The variants provided here can either be in the form of a vcf or merely a list of variants with an appropriate first-line header specifying which column contains the chromosome number, the start position, the reference call and the alternate call. Only tab or comma seperated files of this type are accepted.

-a or --all-variants

REQUIRED. The full path to the vcf file containing all variants. This is usually set to the vcf file generated at the end of the variant calling and filtering pipeline. To only phase selected variants here use a combination of the genomic intervals and filtered variants arguments (-f and -g)

-p or --patient

REQUIRED. The id used to refer to the patient of interest in the vcf and ped files. Must be internally consistent throughout all files.

-o or --output

REQUIRED. The path to the filename where the output should be written. If the file already exists, it will be deleted on program start. Output will be tab seperated format.

-g or --gene-regions

The full path to the bed file containing the genomic regions of interest that should be phased. Only variants within these regions will be phased. Must be in standard bed format.

-r or --reads

