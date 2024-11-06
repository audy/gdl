# `gdl` - Genome Download

[![Rust](https://github.com/audy/gdl/actions/workflows/rust.yml/badge.svg)](https://github.com/audy/gdl/actions/workflows/rust.yml)

A fast, command-line tool for downloading genome assemblies from NCBI. `gdl`
aims to be fast and easy to use.

## Features

- Taxonomy-aware: Fetch all genomes for any node on the NCBI taxonomy tree
  (E.g., fetch all assemblies under the family Enterobacteriaceae)
- Fast: Runs in parallel and is X to Y times faster than other tools

## Installation

### From Source

```sh
git clone ...
cargo build --release
(sudo) cp target/release/gdl /usr/local/bin/
```

## Usage

```
Usage: gdl [OPTIONS] <--tax-id <TAX_ID>|--tax-name <TAX_NAME>>

Options:
      --assembly-summary-path <ASSEMBLY_SUMMARY_PATH>
          path to assembly_summary.txt [default: assembly_summary_refseq.txt]
      --taxdump-path <TAXDUMP_PATH>
          path to extracted taxdump.tar.gz [default: taxdump]
      --dry-run
          do not actually download anything
      --no-cache
          re-fetch assembly_summary.txt and taxdump
      --parallel <PARALLEL>
          [default: 1]
      --format <FORMAT>
          [default: fna] [possible values: fna, faa, gbk, gff]
      --out-dir <OUT_DIR>
          output directory, default=pwd
      --tax-id <TAX_ID>
          tax_id to download assemblies for (includes descendants unless --no-children is enabled)
      --no-children
          do not include child tax IDs of --tax-id (only download assemblies that have the same tax ID as provided by --tax-id)
      --tax-name <TAX_NAME>
          tax_name to download assemblies for (includes descendants unless --no-children is enabled)
      --assembly-level <ASSEMBLY_LEVEL>
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

3. `--cache-dir` - where to store the NCBI tax dump and `assembly_summary.txt`
5. `--no-clobber` - do not overwrite existing files
6. `--verify` - download and check MD5SUM files

- Atomic downloads (don't save partially downloaded files if interrupted)
