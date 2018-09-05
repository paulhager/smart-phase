# SmartPhase

SmartPhase is a phasing tool tailored for clinical use in genetic diagnosis pipelines. It accurately and efficiently reduces the number of possible compound heterozygous variant pairs being examined around either genetic loci or from a list of preselected variant pairs. To achieve this, SmartPhase is able to incorperate parental genotype information as well as reads generated from either DNA- or RNA-sequencing. Furthermore, it incorporates existing haplotype information and applies logical rules to exclude variant constellations that cannot be disease causing.

For a more thorough explanation of SmartPhase and it's validation, please refer to the following: <LINK TO ARTICLE/PREPRINT>

SmartPhase is offered as an executable jar in addition to its source code. The following serves as a general overview and explanation of all its arguments and recommended values.

## Requirements

- Java 8

## Use Case

Clone the repository to your file system and execute the following command from the root folder:

```
java -jar -f ./Validation/UseCase/sim.CEU.trio.chr19.phased.AGV6UTR.hg19GENCODE_codingExons.recode.vcf.gz -a ./Validation/UseCase/sim.CEU.trio.chr19.phased.AGV6UTR.recode.vcf.gz -p NA12878 -r ./Validation/UseCase/simulated.art.hsxt.150l.100fc.400m.100s.NA12878.chr19.bam -m 60 -d ./Validation/UseCase/CEU.ped -o ./Validation/UseCase/smartPhase.sim.NA12878.chr19.trio.results.tsv -x -v -t
```

## Running SmartPhase

SmartPhase can either be run in explorative or analytic mode. In explorative mode, SmartPhase parses provided genomic regions of interest (typically protein coding regions) provided in a bed file and, per region, phases either all variants or merely those specified in a filtered variants file. In analytic or paired mode, no regions of interest must be given, just specific variant pairs that should be phased to one another and SmartPhase will create the appropriate regions. Analytic mode lends itself more to an analysis of a cohort where heavy pre-filtering has already been done and only specific candidate variant pairs should be analysed whereas explorative mode is more all-purpose and can be run at any point in the analysis pipeline.

### Basic usage:

```
java -jar SmartPhase.jar -g /path/to/genomic/regions/to/be/phased/file/allGeneRegionsCanonical.HG19.GRCh37.bed \
-a /path/to/vcf/containing/all/variants/sample.vcf.gz -p PID12345 \
-r /path/to/bams/containing/reads/DNAseq.reads.bam,/path/to/bams/containing/reads/RNAseq.reads.bam -m 60,255 \
-d /path/to/ped/family.ped -t -o /path/to/desired/output/file/output.tsv
```

### Options:

```
-a or --all-variants
```

**REQUIRED**. The path to the vcf file containing all variants. This is usually set to the vcf file generated at the end of the variant calling and filtering pipeline. To only phase selected variants here use a combination of the genomic intervals and filtered variants arguments (-f and -g)

```
-p or --patient
```

**REQUIRED**. The id used to refer to the patient of interest in the vcf and ped files. Must be internally consistent throughout all files.

```
-o or --output
```

**REQUIRED**. The path to the filename where the output should be written. If the file already exists, it will be deleted on program start. Output will be tab seperated format.

```
-f or --filtered-variants
```

The path to the file containing those variants that should be phased. If running in analytic mode, this file is required. If running in explorative mode, this can be left blank and will be set equal to the all-variants file provided. SmartPhase then assumes that all variants that fall within the genomic regions of interest should be phased. The variants provided here can either be in the form of a vcf or merely a list of variants with an appropriate first-line header specifying which column contains the contig, the start position, the reference call and the alternate call. Only tab or comma seperated files of this type are accepted. A common data source for this type of file would for example be a gemini database.

```
-g or --gene-regions
```

The path to the bed file containing the genomic regions of interest that should be phased. Only variants within these regions will be phased. Must be in standard bed format. If this argument is not set, analytic/paired mode is assumed and genomic regions will be built around the variants pairs specified in the filtered variants (-f) file.

```
-r or --reads
```

A comma seperated list of paths pointing to the BAM files containing the reads to be used for read backed phasing. Can either be reads generated through DNAseq or RNAseq. If this option is used, the mapping quality filter option (-m) must also be used and have the same number of items. The order must be the same also. This option and/or trio backed phasing (-t) must be activated.

```
-m or --mapq
```

A comma seperated list of integers indicating the cut-off mapping quality value for reads in the respective file given in the reads flag (-r). For example, if the arguments passed to -r are file1,file2 and the arguments passed to -m are 255,60 this would indicate that only those reads in file1 with a mapping quality of at least 255 are to be considered and only those reads with a mapping quality of at least 60 in file2 are to be considered (in this case most likely because file1 contains reads generated through RNAseq and file2 through DNAseq and only high-confidence, uniquely mapped reads should be used)

```
-t or --trio
```

A boolean switch indicating that trio phasing should be done. If this is present, pedigree information must also be provided (-d). All variant call information for mother, father and child (patient) must be given in the all variants file (-a)

```
-d or --ped
```

The path to the ped file containing information on the familial structure of the trio wishing to be phased. Must be a valid ped file. May contain other families as well as a long as the patient is also present.

```
-y or --physical-phasing
```

Indicates GATK physical phasing information present in the all variants vcf file (-a) should be used as a last resort to phase when no trio or read-backed evidences were available. If using GATK haplotype caller to call variants, this information should already be present. These variants are most likely the result of local haplotype restructuring and thus were not able to be phased using read-backed phasing. In this case, using the information written by GATK HC during restructuring allows SmartPhase to still ascertain phase.

```
-x or --reject-phase
```

A boolean switch indicating that any phase information already annotated in the vcf, for example from other phasing programs, should be ignored. By default this information is taken as correct and those variants that already phased are not phased again. To phase all variants in a file, regardless of other phasing information present, add this argument.

```
-v or --validation
```

A boolean switch indicating that validation files should be generated and output alongside the regular output file. Extra generated files list all haplotype block lengths, the confidence scores all variants correctly and incorrectly called, connectivity information and mean switch error calculations. Files are created in the same location as the regular output file. This flag requires that the true phase of the patient in all variants is in the all variants file(-a).

```
-h or --help
```

Print a summary of the above explanations.

## Paired Mode Input Specification

When running in analytic mode, the filtered variants file (-f) should consist of pre-selected variant pairs for which the phase should be determined. Only tab and comma seperated files are accepted here. The first column may contain multiple comma seperated patient IDs. Only those lines with the patient ID specified in the patient flag (-p) will be phased. If multiple, comma seperated patient IDs are present, the file must be tab seperated. The second and third columns are the two variants whose phase should be determined. Each variant entry must follow the following pattern `contig-start-reference-alternate`. For example: `chr1-143555-G-A`.

A full entry in such a file could look like this:

```
PID12345,4-722294-G-A,4-722315-T-C
```

or like this:

```
PID12345,PID9734,PID2356	4-722294-G-A	4-722315-T-C
```

# Output

## General Specification

The output generated by SmartPhase is stored in a tab-seperated format and consists of 5 columns. 

A typical line in the output file would look like this:

```
C1orf170-1-910578-917473	1-914333-C-G	1-914521-C-T	2	0.7499185705749017
```

The first column is the label of the region this variant pair fell into in the format: `name-contig-start-stop`. The second column shows the first variant in the pair in standard variant format: `contig-start-reference-alternate`. The third column shows the second variant in the pair, also in standard variant format. The fourth column is the flag representing the phasing information determined during execution (see the next section for an explanation of the flags). And the fifth column is a confidence score indicating how confident SmartPhase is in it's phasing call. If a variant pair was phased using trio information, the confidence is set to 1. If read information was used, the confidence takes into account how much conflicting evidence was found and the total number of reads examined. A more thorough explanation of the confidence system can be found in our paper. 

 
## Flags

| Bit | Cis | Trans | Not phased | Innocuous | Not found
|:-----:|:-----:|:-------:|:------------:|:-----------:|:----------:
| 1 | x | | | |
| 2 | | x | | |
| 4 | | | x | |
| 9 | x | | | x |
| 10 | | x | | x |
| 12 | | | x | x |
| 17 | x | | | | x
| 18 | | x | | | x
| 20 | | | x | | x
| 25 | x | | | x | x
| 26 | | x | | x | x
| 28 | | | x | x | x
