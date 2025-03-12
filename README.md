* Annotate structural variants with pangenome graphs

** Generate test data

`misc/generate-test-data.py` is a simple generator that creates a
reference in FASTA and GFA format, and a VCF file with random
variants. You can specify the number and type of variants.

Example use:
```
python misc/generate-test-data.py --output-dir testdata --gfa --vcf --fasta --num-samples 10 --phased -v 20
```
The output is in the `testdata` directory.

** Map VCF variants to GFA

We use GFA to represent variants and haplotypes. Next, we need to map
the variants in the VCF file to paths in a GFA file. There are two
ways to do this. If the VCF file is fully phased, we can map *all*
variants in a sample to two haplotype-specific paths. Example use:
```
python vcf-to-gfa-converter.py testdata/test.gfa testdata/test.vcf o.gfa --allow-mismatches --haplotype-paths --respect-phasing
```
This will create paths in the GFA file that look something like this:
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
python vcf-to-gfa-converter.py testdata/test.gfa testdata/test.vcf o.gfa --allow-mismatches
```
Output will look like this:
```
P       ALT     S1_VAR18_VAR5+,S2_VAR6_VAR14_VAR11+,S3_VAR9_VAR2_VAR1+,S4_VAR17_VAR20_VAR8+,S5_VAR12+,S6+,S7_VAR16_VAR15_VAR13+,S8_VAR10+,S9_VAR3_VAR7_VAR4+,S10_VAR19+,S11+    *       VN:Z:variant_path     SR:Z:REF
```

