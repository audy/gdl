# `gdl`

[![Rust](https://github.com/audy/gdl/actions/workflows/rust.yml/badge.svg)](https://github.com/audy/gdl/actions/workflows/rust.yml)

Genome Download

A fast, command-line tool for downloading genome assemblies from NCBI. `gdl`
aims to be fast, easy to use, and fast.

## Features

- Taxonomy-aware: can fetch all genomes for any node on the NCBI taxonomy tree
- Efficient: starts downloading immediately. No need to pre-fetch thousands of
  unrelated `MD5SUM` files before downloading the first assembly. Written in
  Rust so very low memory consumption

## Installation

TODO

## Usage

```
Usage: gdl [OPTIONS] <--tax-id <TAX_ID>|--tax-name <TAX_NAME>>

Options:
  -a, --assembly-summary-path <ASSEMBLY_SUMMARY_PATH>
          path to assembly_summary.txt [default: assembly_summary_refseq.txt]
  -t, --taxdump-path <TAXDUMP_PATH>
          path to extracted taxdump.tar.gz [default: taxdump]
  -d, --dry-run
          do not actually download anything
  -n, --no-cache
          re-fetch assembly_summary.txt and taxdump
  -p, --parallel <PARALLEL>
          [default: 1]
  -f, --format <FORMAT>
          [default: fna] [possible values: fna, faa, gbk, gff]
  -t, --tax-id <TAX_ID>
          tax_id to download assemblies for (includes descendants)
  -n, --no-children
          do not include child tax IDs of --tax-id (only download assemblies that have the same tax ID as provided by --tax-id)
  -t, --tax-name <TAX_NAME>
          tax_name to download assemblies for (includes descendants)
  -a, --assembly-level <ASSEMBLY_LEVEL>
          include assemblies that match this assembly level. can be used multiple times by default, all assembly_levels are included
  -h, --help
          Print help
```

### Examples

```sh
# download all E. coli genomes (including descendants of E. coli) in FASTA format
gdl --tax-id 562 --format fasta

# download all E. coli genomes in GFF format
gdl --tax-id 562 --format gff
```

# TODO

0. Benchmarks
1. `--format` - fasta, genbank, gff, ...
2. `--out-dir`
3. `--cache-dir` - where to store the NCBI tax dump and `assembly_summary.txt`
4. `--repository` - either GenBank or RefSeq (default is GenBank)
5. `--no-clobber` - do not overwrite existing files
6. `--verify` - download and check MD5SUM files
