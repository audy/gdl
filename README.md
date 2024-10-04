# `fetch-genomes`

A fast, command-line tool for downloading genome assemblies from NCBI.
`fetch-genomes` aims to be fast, easy to use, and fast.

TODO: gif here...

## Features

- Taxonomy-aware: can fetch all genomes for any node on the NCBI taxonomy tree
- Efficient: starts downloading immediately. No need to pre-fetch thousands of
  unrelated `MD5SUM` files before downloading the first assembly. Written in
  Rust so very low memory consumption

## Installation

### Cargo

`cargo install fetch-genomes`

### Homebrew (macOS)

```sh
brew tap audy/fetch-genomes
brew install fetch-genomes
```

### (Windows)

```sh
...
```

### Conda

[Instructions](https://www.theregister.com/2024/08/08/anaconda_puts_the_squeeze_on/)

## Usage

### Download all _Phocaeicola vulgatus_ assemblies in fasta format

Note: this includes all assemblies that are assigned to child nodes of _P. vulgatus_

```sh
fetch-genomes --tax-id 821 --format fasta
```

#### Same as above but only Complete Genomes

```sh
fetch-genomes --tax-id 821 --format fasta --assembly-level "Complete Genome"
```

#### Use GenBank instead of RefSeq

```sh
fetch-genomes --tax-id 821 --format fasta --assembly-level "Complete Genome" --source GenBank
```

#### Run 4 jobs in parallel

```sh
fetch-genomes --tax-id 821 --format fasta --assembly-level "Complete Genome" --source GenBank --n-workers=4
```

### Download all Complete Viral genomes

```sh
fetch-genomes --tax-id 10239 --format fasta --assembly-level "Complete Genome"
```

### Download all human genomes

```sh
fetch-genomes --tax-id 9606 --format fasta --assembly-level "Complete Genome"
```

### Use a taxonomic name instead of a Tax ID

```sh
fetch-genomes --tax-name "Escherichia coli"
```
