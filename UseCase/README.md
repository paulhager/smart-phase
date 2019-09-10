# Use Case

## Underlying Data

Our use case data is based on the simulated data for NA12878 part of the validation data set of SmartPhase as described in the corresponding [README](https://github.com/paulhager/smart-phase/tree/master/Validation).
We selected three pairs of heterozygous missense variants in the genes [ACTRT2](https://omim.org/entry/608535), [CLDN19](https://omim.org/entry/610036) and [PCSK4](https://omim.org/entry/600487). Of each pair one variant has an [gnomAD](https://gnomad.broadinstitute.org/) allele frequency below 5% and the other below 20%. 

## Results

### Standard output/error

The standard output/error informs comprehensively about the phasing process. 
For the use case, there will be messages mentioning that all chromosomes except for chromosome 1 and 19 are skipped.
Finally, there is a summary of the counts of phased and cleared variant pairs in addition to the runtime.

### Output files

The main output file in the `UseCase` folder is `CEU_UseCase_results.tsv`.
If you run SmartPhase in trio mode (-t) the variant pair in ACTRT2 will be phased as being in trans, the pair in CLDN19 will be phased as being in cis and the pair in PCSK4 will be labelled as innocuous.
Assuming that these three pairs are the top candidates for a patient, the variants ACTRT2 should be further be analyzed, the other two genes can be excluded as being pathogenic.
If you run SmartPhase in read-only mode (without the -t option), the variants in ACTRT2 are phased as being in trans and the pair in CLDN19 is phased as being in cis with a high confidence score.
The variants in PCSK4 are too far apart from each other and can therefore not be phased in read-only mode. 

Besides the tabular output file, all phasing information is written to `CEU_UseCase_sp.vcf.gz`.
Each genotype that is phased as part of a haplotype block has annotated the phased genotype (SPGT) and the corresponding haploblock identifier (SPID).
Note that the first genotype of a haplotype is always given as 0|1 as it is not always possible to indicate which is the paternal and the maternal allele.
