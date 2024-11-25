# `gdl` - Genome Download

[![Rust](https://github.com/audy/gdl/actions/workflows/rust.yml/badge.svg)](https://github.com/audy/gdl/actions/workflows/rust.yml)

A fast, easy-to-use, command-line tool for downloading genome assemblies from
NCBI GenBank (RefSeq).

## Features

- **Taxonomic Filtering** - Download all assemblies belonging to a specific
  kingdom, phylum, class, order, family, genus, species, .... Filter based on
  name OR NCBI Tax ID. No need to separately look up tax IDs.
- **Fast** - Use all available processors + connections are re-used to reduce overhead
- **Multiple Sources** - Fetch assemblies from RefSeq or GenBank
- **Filtering** - Filter based on `assembly_level` (Complete, Contig, Chromosome, ...)

## Examples

```sh
# Download all complete genomes in Lactobacillales (order) in GenBank format
gdl --tax-name "Lactobacillales" --format gbk --source refseq --assembly-level "Complete Genome"

# Download all Betacoronavirus genomes currently in GenBank
gdl --tax-name "Betacoronavirus" --format fna --source genbank --out-dir betacoronaviruses/

# Download all Complete viral assemblies
gdl --tax-name "Viruses" --format fna --source refseq --assembly-level "Complete Genome"
```

## Installation

[Download a pre-compiled binary]()

### From Source

```sh
git clone ...
cargo build --release
(sudo) cp target/release/gdl /usr/local/bin/
```

### From Cargo

```sh
cargo install gdl
```

## Full Usage

```
Usage: gdl [OPTIONS] <--tax-id <TAX_ID>|--tax-name <TAX_NAME>>

Options:
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
      --source <SOURCE>
          where to fetch assemblies from (default is RefSeq) [default: refseq] [possible values: genbank, refseq, none]
      --assembly-summary-path <ASSEMBLY_SUMMARY_PATH>
          path to assembly_summary.txt
      --tax-id <TAX_ID>
          tax_id to download assemblies for (includes descendants unless --no-children is enabled)
      --no-children
          do not include child tax IDs of --tax-id (only download assemblies that have the same tax ID as provided by --tax-id)
      --tax-name <TAX_NAME>
          tax_name to download assemblies for (includes descendants unless --no-children is enabled)
      --assembly-level <ASSEMBLY_LEVEL>
          include assemblies that match this assembly level. By default, all assembly_levels are included
  -h, --help
          Print help
```
