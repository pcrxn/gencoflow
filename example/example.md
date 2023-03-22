# Example

1. The genomic assembly for a pathogenic Escherichia coli isolate (BioSample SAMN02147037) was retrieved from https://zenodo.org/record/5604579, and corresponds to the following assembly methods described in the paper by Raphenya et al. 2022:

- Filtered/GenBank/Raw: Filtered
- Assembler: SKESA

More metadata for this isolate is available at https://zenodo.org/record/6543963 under `data/strain.metadata.txt`. 

2. Prokka (v1.14.6) was used to identify and annotate ORFs within the assembly using the following commands:

```bash
$ cd gencoflow/examples/
$ prokka --cpus 14 --outdir prokka_SAMN02147037 SAMN02147037.fasta
```

All output files except for the GenBank file were removed.

3. ABRicate (v1.0.1) was used to identify antimicrobial resistance genes within the assembly using the following commands:

```bash
$ cd gencoflow/examples/
$ abricate --db card SAMN02147037.fasta > card_SAMN02147037.tsv
Using nucl database card:  2631 sequences -  2022-Nov-25
Processing: SAMN02147037.fasta
Found 45 genes in SAMN02147037.fasta
Tip: found a bug in abricate? Post it at https://github.com/tseemann/abricate/issues.
Done.
```

4. Columns in card_SAMN02147037.tsv were renamed to the following: 'seq_id', 'contig_id', 'start', 'end', 'strand', 'arg', 'cov', 'cov_map', 'gaps', 'percent_cov', 'pident', 'database', 'accession', 'product', 'resistance'.

## References

Raphenya AR, Robertson J, Jamin C, de Oliveira Martins L, Maguire F, McArthur AG, Hays JP. 2022. Datasets for benchmarking antimicrobial resistance genes in bacterial metagenomic and whole genome sequencing. Sci Data 9:341. doi: 10.1038/s41597-022-01463-7
