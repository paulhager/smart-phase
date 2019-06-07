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

The main output file is the `UseCase` folder is `smartPhase.sim.NA12878.chr19.muc16.trio.results.tsv`.

| Bitflag | Interpretation | Count |
|:-------:|:--------------:|:-----:|
| 1 | Cis | 1,035 |
| 2 | Trans | 46 |
| 9 | Innocuous & Cis | 542 |
| 10 | Innocuous & Trans | 1,458 |
| 17 | Cis & Not found | 5,409 |
| 18 | Trans & Not found | 296 |
| 25 | Innocuous & Cis & Not found | 3,790 |
| 26 | Innocuous & Trans & Not found | 7,725 |
| 28 | Innocuous & Not Phased & Not found | 202 |
