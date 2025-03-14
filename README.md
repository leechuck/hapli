# hapli - Haplotype and Genotype Functional Annotation Tool

hapli provides functional annotations to haplotypes and genotypes using GFA files as a basis for annotation.

## Installation

### Quick Install

```bash
pip install hapli
```

### Development Install

Clone the repository and install in development mode:

```bash
git clone https://github.com/leechuck/hapli.git
cd hapli
pip install -e .
```

### Install with Optional Dependencies

```bash
# Install with RDF support
pip install -e .[rdf]

# Install with sequence alignment support
pip install -e .[align]

# Install with development dependencies
pip install -e .[dev]

# Install all optional dependencies
pip install -e .[all]
```

### From Requirements File

Install dependencies using the requirements file:

```bash
pip install -r requirements.txt
```

### Requirements

- Python 3.8 or higher
- BioPython (for sequence analysis)
- NumPy (for numerical operations)
- PyYAML (for configuration)
- Pandas (for data manipulation)

#### Optional Dependencies
- RDFLib (for RDF output)
- NetworkX (for graph operations)
- Parasail and Edlib (for efficient sequence alignment)

## License

This project is licensed under the GNU General Public License v3 (GPLv3).

## Usage

hapli has several commands for different tasks:

### Analyze Variant Effects

```bash
# Basic analysis
python hapli.py analyze input.gfa annotations.gff3

# Generate RDF output
python hapli.py analyze input.gfa annotations.gff3 --format rdf --output output.ttl

# Generate sample-specific reports
python hapli.py analyze input.gfa annotations.gff3 --sample-reports --output-prefix sample_reports

# Use sequence alignment for more accurate variant effect prediction
python hapli.py analyze input.gfa annotations.gff3 --use-alignment

# Generate a consolidated RDF report for all samples
python hapli.py analyze input.gfa annotations.gff3 --format rdf --consolidated --output consolidated.ttl
```

### Generate Test Data

```bash
# Generate all test data types
python hapli.py generate-test-data --output-dir testdata --gfa --gff3 --vcf --fasta

# Customize test data generation
python hapli.py generate-test-data --output-dir testdata --length 5000 --segments 10 --variants 20 --genes 15 --samples 5
```

### Convert VCF to GFA

```bash
# Basic conversion
python hapli.py vcf-to-gfa input.gfa input.vcf output.gfa

# Allow reference mismatches
python hapli.py vcf-to-gfa input.gfa input.vcf output.gfa --allow-mismatches

# Create haplotype-specific paths with phasing information
python hapli.py vcf-to-gfa input.gfa input.vcf output.gfa --haplotype-paths --respect-phasing
```

### Generate GFF3 from GFA

```bash
# Basic conversion
python hapli.py gfa-to-gff3 input.gfa output.gff3

# Specify number of genes and exons
python hapli.py gfa-to-gff3 input.gfa output.gff3 --genes 25 --exons 3
```

## Example Output

The analyze command generates detailed reports about variant effects on genomic features:

```
# Variant Effect Report for Sample: SAMPLE1
#===============================================================================

## Variant Summary
--------------------------------------------------------------------------------
No variants identified for this sample.


## Feature Effect Summary for Sample: SAMPLE1
--------------------------------------------------------------------------------
Total Features Analyzed: 40
Features Affected by Variants: 30

Zygosity Summary:
  Homozygous Effects: 14
  Heterozygous Effects: 16

Affected Features by Type:
  CDS: 15/15 (100.0%)
  exon: 15/15 (100.0%)

Variant Effects by Type:
  splicing_affected: 15
  amino_acid_change: 15
  in_frame_change: 15
  start_codon_disruption: 5
  premature_stop_codon: 5


## Detailed Feature Effects for Sample: SAMPLE1

### CDS Features
--------------------------------------------------------------------------------

Feature: cds1.1 (cds1.1)
Location: 291-364 (+)
Zygosity: HOMOZYGOUS
Effects:
  - amino_acid_change
  - in_frame_change
  - start_codon_disruption

Amino Acid Changes:
  - H1I
  - T2H
  - *3T
  - V4*
  - T5V
  ...
```

## Features

- Parse GFA files with variant information
- Parse GFF3 files with genomic feature annotations
- Analyze variant effects on genomic features
- Generate detailed reports in text or RDF format
- Support for phased haplotypes and compound heterozygous variants
- Sequence alignment-based variant effect prediction
- Generate test data for development and testing
- Convert between different genomic data formats

## License

MIT
