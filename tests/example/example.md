# Example

1. The genomic assemblies for several pathogenic isolates (Escherichia coli, BioSample SAMN02147037; Acinetobacter baumannii, BioSample SAMN10163228; Salmonella enterica, BioSample SAMN05263515) were retrieved from https://zenodo.org/record/5604579, and correspond to the following assembly methods described in the paper by Raphenya et al. 2022:

- Filtered/GenBank/Raw: Filtered
- Assembler: SKESA

More metadata for these isolates are available at https://zenodo.org/record/6543963 under `data/strain.metadata.txt`. 

The assemblies were saved as .fastq.gz files within `examples/`.

2. Prokka (v1.14.6) was used to identify and annotate ORFs within the assemblies using default parameters.

All output files except for the GenBank files were removed.

3. ABRicate (v1.0.1) was used to identify antimicrobial resistance genes within the assemblies using the pre-loaded CARD database. For multiple samples simultaneously:

```bash
$ cd examples/
$ abricate --db card *.fasta > multi-sample/card-results.tsv
```

The '.fasta' at the end of each sample ID in the '#FILE' column of `card-results.tsv` was removed:

```bash
sed -i 's/.fasta//' multi-sample/card-results.tsv
```

## References

Raphenya AR, Robertson J, Jamin C, de Oliveira Martins L, Maguire F, McArthur AG, Hays JP. 2022. Datasets for benchmarking antimicrobial resistance genes in bacterial metagenomic and whole genome sequencing. Sci Data 9:341. doi: 10.1038/s41597-022-01463-7
