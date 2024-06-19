# Nextflow Pipeline for Metabolome-Wide Colocalization Analysis

## Introduction

This pipeline is a variation of the `nf-scallop-coloc` pipeline available
[here](https://github.com/comp-med/nf-scallop-coloc). Instead of running
colocalization analysis of pQTL data on CVD outcomes, we use the UK Biobank NMR
data mQTL results. A main difference is that this pipeline requires an existing
LD matrix for each region and does not compute one manually.

## Input Data

## Requirements

* An `apptainer` image used for all processes that execute `R` code
* A local library containing the `R` packages used in the analysis:
    * `data.table`
    * `glue`
    * `biomaRt`
    * `susieR`
    * `coloc`
    * `doMC`
    * `Rfast`
* Various UK Biobank related files outlined below

## External Input Files

The pipeline requires various external inputs with specific specifically named
files to function. These directories need to be specified in the `params` block
in the file `nextflow.config` (or the configuration for the used profile in
`config/`).

### UK Biobank NMR mQTL Finemapping Results 

The main input from the mQTL analysis are the finemapped association results
for each region. To process all input files correctly, a master table needs to
specified with the parameter `params.nmr_finemap_master_table`. It contains
information about each metabolite, region, and where the results are located.

```
index  phenotype           region_name  chrom  region_start  region_end
1      Ala                 1_0          1      603426        2638126
2      Val                 1_0          1      603426        2638126
3      Total_BCAA          1_0          1      603426        2638126
4      Citrate             1_0          1      603426        2638126
5      ApoA1               1_0          1      603426        2638126
[...]
```

The directory containing the files for each metabolite is to be specified using
the parameter `params.nmr_credible_sets_directory`. Results for each
combination and metabolite are saved in sequentially numbered sub-directories.
The content of the directory looks as follows:

```bash
.
├── 1
│   ├── finemapped_results.txt.gz
│   ├── flt_dosages.bgen
│   ├── snps.dosage
│   └── variants_dosages.txt
├── 10
│   ├── finemapped_results.txt.gz
│   ├── flt_dosages.bgen
│   ├── snps.dosage
│   └── variants_dosages.txt
[...]
```

The pre-computed LD matrix is expected to be located in a directory containing
region-specific data. The main directory is specified using the parameter
`params.nmr_ld_directory`. The the root directory structure is as follows.

```bash
.
├── 10_0
│   ├── ldmat.ld
│   ├── mfile.z
│   ├── snpdata.z
│   └── snplist.txt
├── 10_1
│   ├── ldmat.ld
│   ├── mfile.z
│   ├── snpdata.z
│   └── snplist.txt
[...]
```

The files themselves are custom and all analyses require that the column names
adhere to the expected format.

### Recombination Map

For plotting colocalization results, files containing recombination rates in cM
are required for each chromosome. The directory containing these files is to be
specified in `params.recombination_rate_map_dir`. The files in the directory
must to be named as seen below:

```bash
.
├── genetic_map_GRCh37_chr1.txt
├── genetic_map_GRCh37_chr2.txt
├── genetic_map_GRCh37_chr3.txt
├── genetic_map_GRCh37_chr4.txt
├── genetic_map_GRCh37_chr5.txt
├── genetic_map_GRCh37_chr6.txt
├── genetic_map_GRCh37_chr7.txt
├── genetic_map_GRCh37_chr8.txt
├── genetic_map_GRCh37_chr9.txt
[...]
└── genetic_map_GRCh37_chrX.txt
```

The files must adhere to the following structure:

```bash
Chromosome	Position(bp)	Rate(cM/Mb)	Map(cM)
chr1	55550	2.981822	0.000000
chr1	82571	2.082414	0.080572
chr1	88169	2.081358	0.092229
chr1	254996	3.354927	0.439456
chr1	564598	2.887498	1.478148
chr1	564621	2.885864	1.478214
chr1	565433	2.883892	1.480558
chr1	568322	2.887570	1.488889
chr1	568527	2.895420	1.489481
[...]
```

### Outcome directory

The outcome data must be pre-processed using
[`MungeSumstats`](https://neurogenomics.github.io/MungeSumstats/). The
directory containing the statistics is specified in the parameter
`params.outcome_sumstat_directory`. The files are organized in this directory
as follows:

```bash
.
├── phenotype_1
│   ├── GRCh37
│   │   └── phenotype_1_grch37.tsv.bgz
│   └── GRCh38
│       └── phenotype_1_grch38.tsv.bgz
└── phenotype_2
    ├── GRCh37
    │   └── phenotype_2_grch37.tsv.bgz
    └── GRCh38
        └── phenotype_2_grch38.tsv.bgz
```

The data dictionary `params.outcome_data_dictionary` contains meta data for
each phenotype. Specifically, the sample size is necessary for the
colocalization analysis.

## Get Started

* Add parameters in `nextflow.config`
* Add container information in `config/cluster.config`

Run the pipeline with:

```bash
# To run the pipeline using SLURM and apptainer containers from the local HPC
nextflow run main.nf -resume -profile cluster

# To get process output add `-process.debug`
nextflow run main.nf -resume -process.debug -profile cluster
```

## Trouble Shooting

When jobs fail, you can go to the work directory with the failed job and use
`srun` or `sbatch` to manually submit the job and check the output.

```bash
cd work/<PATH/TO/FAILED/JOB>

# Manually submit the job
srun --mem '64G' --ntasks=1 --cpus-per-task=8 .command.run
```
