# Use Case

## Underlying Data

Our use case data is based on the simulated for the validation of SmartPhase as described in the corresponding [README](https://github.com/paulhager/smart-phase/tree/master/Validation).
We extracted all variants identified in the CEU trio in the gene MUC16 to create the VCF file and extracted all simulated reads mapping to MUC16 to create the BAM file.
In total, there are 203 heterozygous variants in the child which results in 20,503 potential compound heterozygous variant pairs.

## Results

### Standard output/error

The standard output/error informs comprehensively about the phasing process. 
For the use case, there will be quite some messages mentioning that all chromosomes except of chromosome 19 are skipped.
Then, there are error messages that some variants in the VCF cannot be found in the aligned reads in the BAM file.
Finally, there is a summary of the counts of phased and cleared variant pairs and the runtime is given.

### Output files

The main output file in the `UseCase` folder is `smartPhase.sim.NA12878.chr19.muc16.trio.results.tsv`.
The table below summarizes the resulting bitflags for the 20,503 variant pairs.
All phasing predictions resulted from trio-phasing with a confidence score of 1 and are in concordance with the simulated haplotypes.

| Bitflag | Interpretation | Count |
|--------:|----------------|------:|
| 1 | Cis | 1,035 |
| 2 | Trans | 46 |
| 9 | Innocuous & Cis | 542 |
| 10 | Innocuous & Trans | 1,458 |
| 17 | Cis & Not found | 5,409 |
| 18 | Trans & Not found | 296 |
| 25 | Innocuous & Cis & Not found | 3,790 |
| 26 | Innocuous & Trans & Not found | 7,725 |
| 28 | Innocuous & Not Phased & Not found | 202 |

Under the assumption that the child suffers from a disease and both parents are healthy, there are only 342 compound heterozygous variant that are potentially disease causing.
All other pair are either phased as being located on the same allele or labeled as innocuous variant constellations.

Besides the tabular output file all phasing information is written to `sim.CEU.trio.chr19.phased.muc16_sp.vcf.gz`.
Each genotype that is phased as part of a haplotype block has annotated the phased genotype (SPGT) and the corresponding haploblock identifier (SPID).
For the use case, the genotypes of 202 variants are phased as part of the haplotype block 19\_8962315.
Note that the first genotype of a haplotype is always given as 0|1 as it is not always possible to indicate which is the paternal and the maternal allele.
In our case, one variant cannot be phased because the corresponding genotypes in the both parents are also heterozygous.
All 202 pairs that can be formed with this variant are labeled as innocuous (flag 28 in the table above).


