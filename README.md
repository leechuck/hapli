* Annotate structural variants with pangenome graphs

** Generate test data

`misc/generate-test-data.py` is a simple generator that creates a
reference in FASTA and GFA format, and a VCF file with random
variants. You can specify the number and type of variants.

Example use:
```
python misc/generate-test-data.py --output-dir testdata --gfa --vcf --fasta -v 20
```
