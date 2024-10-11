# `gdl`

Genome Download

A fast, command-line tool for downloading genome assemblies from NCBI. `gdl`
aims to be fast, easy to use, and fast.

# TODOs

1. `--format` - fasta, genbank, gff, ...
2. `--out-dir`
3. `--cache-dir` - where to store the NCBI tax dump and `assembly_summary.txt`
4. `--repository` - either GenBank or RefSeq (default is RefSeq)
5. `--no-clobber` - skip existing files
6. `--verify` - download and check MD5SUM files
7. `--dry-run` - do not actually download anything
8. A single, global progress bar with an ETA

## Features

- Taxonomy-aware: can fetch all genomes for any node on the NCBI taxonomy tree
- Efficient: starts downloading immediately. No need to pre-fetch thousands of
  unrelated `MD5SUM` files before downloading the first assembly. Written in
  Rust so very low memory consumption

## Installation

TODO

## Usage

### Examples

```sh
# download all E. coli genomes in FASTA format
gdl --tax-id 562 --format fasta

# download all E. coli genomes in GFF format
gdl --tax-id 562 --format gff
```

### Advanced Filtering

```sh
gdl \
    --tax-id 512 \
    --tax-id 666 \
    --include assembly_summary 'Complete*' \ # glob match any string field
    --include organism_name "*foo*" \        # case insensitive, glob match
    --include organism_name "/.*foo.*/" \    # case sensitive as defined by regex
    --include organism_name "/.*foo.*/i" \   # case insensitive
    --include genome_size '<50000000' \      # integer filter
    --include seq_rel_date '>2024' \         # not sure about this one yet...
    --exclude organism_name '*foo*' \        # all filters are combined with OR or OR NOT (in the case of exclude)
    --limit 3 \
    --limit-per-species 3 \
    --sort-by assembly_level \ # orders by Complete Genome > Scaffold > Contig > Other
    asdf
```

## Field Names in `assembly_summary.txt:

```
assembly_accession
bioproject
biosample
wgs_master
refseq_category
taxid
species_taxid
organism_name
infraspecific_name
isolate
version_status
assembly_level
release_type
genome_rep
seq_rel_date
asm_name
asm_submitter
gbrs_paired_asm
paired_asm_comp
ftp_path
excluded_from_refseq
relation_to_type_material
asm_not_live_date
assembly_type
group
genome_size
genome_size_ungapped
gc_percent
replicon_count
scaffold_count
contig_count
annotation_provider
annotation_name
annotation_date
total_gene_count
protein_coding_gene_count
non_coding_gene_count
pubmed_id
```
