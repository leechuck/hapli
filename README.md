# Annotate structural variants with pangenome graphs

## Generate test data

`misc/generate-test-data.py` is a simple generator that creates a
reference in FASTA and GFA format, and a VCF file with random
variants. You can specify the number and type of variants.

Example use:
```
python misc/generate-test-data.py --output-dir testdata --gfa --vcf --fasta --num-samples 10 --phased -v 20
```
The output is in the `testdata` directory.

We also create some example genome annotations in GFF3 format. We will
use this to determine the functional impact of the (structural)
variants we identify. The script takes the reference path in a GFA
file as input and makes up some genomic features.

Example usage:
```
python misc/gfa-to-gff3.py --genes 25 --exons 2 testdata/test.gfa testdata/test.gff3
```
This will create a GFF3 file which looks like this:
```
##gff-version 3
##sequence-region 1 1 1000000
1       gfa_to_gff3     promoter        1       35      .       +       .       ID=promoter_gene3;Name=promoter_gene3;Associated_gene=gene3
1       gfa_to_gff3     promoter        1       55      .       +       .       ID=promoter_gene4;Name=promoter_gene4;Associated_gene=gene4
1       gfa_to_gff3     exon    1       83      .       +       .       ID=exon2.1;Parent=mRNA2;Name=gene_2_exon1
1       gfa_to_gff3     CDS     1       83      .       +       .       ID=cds2.1;Parent=mRNA2;Name=gene_2_cds1
1       gfa_to_gff3     exon    1       123     .       -       .       ID=exon1.1;Parent=mRNA1;Name=gene_1_exon1
1       gfa_to_gff3     CDS     1       123     .       -       .       ID=cds1.1;Parent=mRNA1;Name=gene_1_cds1
1       gfa_to_gff3     promoter        1       135     .       +       .       ID=promoter_gene6;Name=promoter_gene6;Associated_gene=gene6
1       gfa_to_gff3     gene    1       262     .       -       .       ID=gene1;Name=gene_1;locus_tag=LOCUS0001
1       gfa_to_gff3     mRNA    1       262     .       -       .       ID=mRNA1;Parent=gene1;Name=gene_1_transcript
1       gfa_to_gff3     gene    1       434     .       +       .       ID=gene2;Name=gene_2;locus_tag=LOCUS0002
1       gfa_to_gff3     mRNA    1       434     .       +       .       ID=mRNA2;Parent=gene2;Name=gene_2_transcript
1       gfa_to_gff3     terminator      36      89      .       -       .       ID=terminator_gene5;Name=terminator_gene5;Associated_gene=gene5
1       gfa_to_gff3     exon    36      256     .       +       .       ID=exon3.1;Parent=mRNA3;Name=gene_3_exon1
1       gfa_to_gff3     CDS     36      256     .       +       .       ID=cds3.1;Parent=mRNA3;Name=gene_3_cds1
1       gfa_to_gff3     gene    36      368     .       +       .       ID=gene3;Name=gene_3;locus_tag=LOCUS0003
1       gfa_to_gff3     mRNA    36      368     .       +       .       ID=mRNA3;Parent=gene3;Name=gene_3_transcript
1       gfa_to_gff3     exon    56      108     .       +       .       ID=exon4.1;Parent=mRNA4;Name=gene_4_exon1
1       gfa_to_gff3     CDS     56      108     .       +       .       ID=cds4.1;Parent=mRNA4;Name=gene_4_cds1
1       gfa_to_gff3     gene    56      325     .       +       .       ID=gene4;Name=gene_4;locus_tag=LOCUS0004
1       gfa_to_gff3     mRNA    56      325     .       +       .       ID=mRNA4;Parent=gene4;Name=gene_4_transcript
...
```

## Map VCF variants to GFA

We use GFA to represent variants and haplotypes. Next, we need to map
the variants in the VCF file to paths in a GFA file. There are two
ways to do this. If the VCF file is fully phased, we can map *all*
variants in a sample to two haplotype-specific paths. Example use:
```
python vcf-to-gfa-converter.py testdata/test.gfa testdata/test.vcf out.gfa --allow-mismatches --haplotype-paths --respect-phasing
```
This will create paths in the `out.gfa` file that look something like this:
```
P       SAMPLE_Sample1_HAP1_1   S1_VAR18+,S2+,S3+,S4_VAR17_VAR20_VAR8+,S5+,S6+,S7_VAR15_VAR13+,S8+,S9_VAR7+,S10+,S11+   *       SM:Z:Sample1    PH:Z:haplotype_1        VA:Z:VAR7;VAR15;VAR13;VAR17;VAR20;VAR8;VAR18  SR:Z:REF
P       SAMPLE_Sample1_HAP2_1   S1_VAR18_VAR5+,S2_VAR14_VAR11+,S3_VAR9_VAR2_VAR1+,S4_VAR17_VAR20_VAR8+,S5+,S6+,S7_VAR15_VAR13+,S8_VAR10+,S9_VAR3_VAR7_VAR4+,S10_VAR19+,S11+     *       SM:Z:Sample1  PH:Z:haplotype_2        VA:Z:VAR19;VAR3;VAR7;VAR4;VAR10;VAR15;VAR13;VAR17;VAR20;VAR8;VAR9;VAR2;VAR1;VAR14;VAR11;VAR18;VAR5      SR:Z:REF
```

However, many VCF files are not fully phased, and we cannot create one
path per haplotype; if we have unphased VCF files, we create two paths
per variant per sample:
```
python vcf-to-gfa-converter.py testdata/test.gfa testdata/test.vcf o.gfa --allow-mismatches --haplotype-paths
```
The output will look something like this:
```
P       SAMPLE_Sample1_VAR_VAR19        S1+,S2+,S3+,S4+,S5+,S6+,S7+,S8+,S9+,S10_VAR19+,S11+     *       SM:Z:Sample1    PH:Z:individual VA:Z:VAR19      SR:Z:REF
P       SAMPLE_Sample1_VAR_VAR3 S1+,S2+,S3+,S4+,S5+,S6+,S7+,S8+,S9_VAR3+,S10+,S11+      *       SM:Z:Sample1    PH:Z:individual VA:Z:VAR3       SR:Z:REF
P       SAMPLE_Sample1_VAR_VAR7 S1+,S2+,S3+,S4+,S5+,S6+,S7+,S8+,S9_VAR7+,S10+,S11+      *       SM:Z:Sample1    PH:Z:individual VA:Z:VAR7       SR:Z:REF
P       SAMPLE_Sample1_VAR_VAR4 S1+,S2+,S3+,S4+,S5+,S6+,S7+,S8+,S9_VAR4+,S10+,S11+      *       SM:Z:Sample1    PH:Z:individual VA:Z:VAR4       SR:Z:REF
P       SAMPLE_Sample1_VAR_VAR10        S1+,S2+,S3+,S4+,S5+,S6+,S7+,S8_VAR10+,S9+,S10+,S11+     *       SM:Z:Sample1    PH:Z:individual VA:Z:VAR10      SR:Z:REF
```
This is closest to the representation of variants in traditional
unphased VCF files, one variant per line.

We can also merge all variants into a single path, e.g., if we want to
do some population-level analysis:
```
python vcf-to-gfa-converter.py testdata/test.gfa testdata/test.vcf out.gfa --allow-mismatches
```
Output will look like this:
```
P       ALT     S1_VAR18_VAR5+,S2_VAR6_VAR14_VAR11+,S3_VAR9_VAR2_VAR1+,S4_VAR17_VAR20_VAR8+,S5_VAR12+,S6+,S7_VAR16_VAR15_VAR13+,S8_VAR10+,S9_VAR3_VAR7_VAR4+,S10_VAR19+,S11+    *       VN:Z:variant_path     SR:Z:REF
```

## Variant effect annotation

Next, we can get "variant" effects by overlaying the path information
generated from the variants with a GFF3 file that contains genomic
features. However, in contrast to methods like VEP, this will not
generate reports for single "variants" but rather haplotype- and
sample-level reports, essentially combining multiple haplotypes and
determine which genomic features are affected by it.

Example use:
```
python variant-effect-report.py out.gfa testdata/test.gff3 --sample-reports --consolidated
```
This will generate test output for all samples in the GFA file, and
will look something like this:
```
# Variant Effect Report for Sample: Sample1
#===============================================================================

## Variant Summary
--------------------------------------------------------------------------------
Variant: VAR13 (INS)
Position: 60-169
Length Change: +9 bp
Affected Segments: S60_VAR13

Variant: VAR1 (DEL)
Position: 59-149
Length Change: -10 bp
Affected Segments: S59_VAR1

Variant: VAR20 (SNP)
Position: 53-159
Length Change: 0 bp (substitution)
Affected Segments: S53_VAR20

Variant: VAR16 (SNP)
Position: 25-115
Length Change: 0 bp (substitution)
Affected Segments: S25_VAR16

Variant: VAR3 (SNP)
Position: 24-121
Length Change: 0 bp (substitution)
Affected Segments: S24_VAR3

Variant: VAR12 (DEL)
Position: 7-106
Length Change: -70 bp
Affected Segments: S7_VAR12, S7_VAR12_VAR18_VAR15, S7_VAR12_VAR18

Variant: VAR19 (INS)
Position: 85-222
Length Change: +38 bp
Affected Segments: S85_VAR19

Variant: VAR4 (INS)
Position: 73-188
Length Change: +8 bp
Affected Segments: S73_VAR4_VAR14

...
```

Of course you can also get the output in any RDF format you like by
appending `--format rdf` to the command (and optionally ``--rdf-format
{turtle,n3,xml,json-ld,ntriples}`).

The "variant types" are identified based on the sequence of nodes in
the path in the GFA file, essentially using positional
information. This is not what we would want because multiple nodes may
change position, and a positional approach will not work
correctly. Instead, we can align the features to the generated
haplotype to identify what the exact changes are:
```
python variant-effect-report.py out.gfa testdata/test.gff3 --sample-reports  --consolidated --use-alignment
```
This will generate a similar report, but identifies the feature
changes by alignments, not just by the sequence of nodes in the path.

