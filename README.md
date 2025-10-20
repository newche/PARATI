# PARATI: PARental Allele Transmission Inference
PARATI infers maternal and paternal transmitted and non‑transmitted alleles from trio genotype data, enabling SNP‑specific analyses of genetic nurture and transgenerational effects.


# Installation

Requirements:
OS: Linux or macOS
Shell: bash
R: ≥ 4.1 (tested with current R release)
R packages: dplyr (1.1.4), data.table (1.16.0), stringr (1.5.1), optparse (1.7.5), openxlsx (4.2.7.1), vcfR (1.15.0)
Optional: PLINK 1.9 on PATH if you want PLINK‑format output

Get the code & install R packages
```
# clone
git clone https://github.com/newche/PARATI.git
cd PARATI
```


# Inputs

PARATI needs two files:

**1. Trio genotype VCF**

Format: .vcf.gz with standard columns CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, followed by sample columns (individual IDs).
Contains biological father, mother, and child for each trio.
Autosomes & biallelic SNPs only.
The VCF does not need to be chromosome‑specific; PARATI works per chromosome internally.
Important: Sample IDs in the VCF must exactly match the IndividualID entries in the family index (case‑sensitive). To avoid parsing issues, do not use underscores _ in sample IDs.
Quality control (recommended before running PARATI): Filter by call rate and MAF; remove Mendelian inconsistencies within trios; restrict to autosomal biallelic SNPs.


Example: 
The following is an example for the Trio genotype VCF input, part of the simulated testing set sim_trios_dash_5k_22chr_single.vcf.gz. 
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1-M	1-P	1-B	2-M	2-P	2-B
1	100000	rs1	G	T	.	PASS	.	GT	0/0	0/1	0/1	0/0	0/1	0/1
1	101000	rs2	A	C	.	PASS	.	GT	0/0	0/1	0/0	0/0	0/0	0/0
1	102000	rs3	A	T	.	PASS	.	GT	0/0	0/0	0/0	1/1	1/1	0/1
1	103000	rs4	A	G	.	PASS	.	GT	0/1	0/1	0/1	0/1	0/0	./.
1	104000	rs5	A	G	.	PASS	.	GT	0/1	0/0	0/1	0/0	0/1	0/0
1	105000	rs6	T	C	.	PASS	.	GT	0/0	1/1	0/1	1/1	0/1	0/1
1	106000	rs7	T	C	.	PASS	.	GT	0/0	0/0	0/0	0/1	0/0	0/1
1	107000	rs8	A	C	.	PASS	.	GT	0/0	0/0	0/0	./.	0/1	0/1
...
```

**2. Family index**

Format: .xlsx with three required columns:
FamilyIndex — family ID, can be integers or characters (e.g., FAM001)
IndividualID — must match a VCF sample column exactly
Role — one of F (father), M (mother), C (child)

Example:
```
FamilyIndex	IndividualID	Role
1	1-M	M
1	1-P	F
1	1-B	C
2	2-M	M
2	2-P	F
2	2-B	C
3	3-M	M
3	3-P	F
3	3-B	C
...
```


# Quick start

**Minimal run** (VCF outputs):

```
Rscript path/to/PARATI/PARATI_script.R \
  --geno path/to/trios.vcf.gz \
  --family path/to/family_index.xlsx \
  --chr 22 \
  --out path/to/results/
```

**With haplotype window tuning, save intermediates, and PLINK outputs**:

```
Rscript path/to/PARATI/PARATI_script.R \
  --geno path/to/trios.vcf.gz \
  --family path/to/family_index.xlsx \
  --chr 1 \
  --out path/to/results/ \
  --haplen 500 \
  --savetemp T \
  --makebed T
```

--chr should match your VCF chromosome naming convention (e.g., 1 vs chr1).
--haplen is the length in kilobases of haplotype to infer heterozygotes (default ±500 kb).


**Command‑line options:**

Required:
--geno : Trio genotype file in .vcf.gz with full path
--family : Family index .xlsx with full path
--chr : chromosome to process (e.g., 1, 2, …, 22)
--out : output directory (created if missing)

Optional:
--haplen <kb> : haplotype half‑window for triple‑heterozygote inference (default 500 → ±500 kb)
--savetemp T|F : keep intermediate files, including a summary table for haplotype calls (default F)
--makebed T|F : also export PLINK .bed/.bim/.fam (default F)



# Outputs

PARATI writes two genotype datasets per run:
1. Parental non‑transmitted genotypes
2. Parental transmitted genotypes

Each dataset:
* Contains the requested chromosome only.
* Preserves SNP metadata (VCF) and uses the parental IndividualIDs from your family index.
* Is saved in VCF format; if --makebed T, also saved as PLINK files.

Note on PLINK output: Because the output encodes only the transmitted or non‑transmitted allele per SNP, alleles are represented as homozygous in standard PLINK format. When computing PRS from these files, use standardized scores (e.g., z‑scores) rather than raw allele counts.

# Testing Example:

```
git clone https://github.com/newche/PARATI.git
cd PARATI

module load plink/1.9

Rscript PARATI_script.R --geno ./test/Toy_TrioGenotype.vcf.gz --family ./test/Toy_FamilyIndexTable.xlsx --out ./ --chr 22 --savetemp T --makebed T

```





