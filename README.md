# scampiman
A pipeline to align and summarize tiled amplicon coverage (of a virus, probably) from raw sequences

**Requires: Reads, Reference Genome(s), Primer `.bed` File**

**Produces: Samtools Ampliconstats File, Table of Amplicon Coverage `.tsv`**

 - Input formats
    `.fastq` or `.bam`
 - Input logic
    `file(s)` or `directory` (with files)
 - Seq tech
    `illumina short read` or `ONT`

1) Align reads to reference (`dorado aligner` or `minimap2`)
2) samtools: view, sort, ampliconclip, ampliconstats
3) Parse ampliconstats output into table, output `.tsv`


# Install

clone this repo

## Make conda environment

`conda env create -f scampiman/environment/scampiman.yml`

once complete:

`conda activate scampiman`

## pip install scampiman

Using `pip` to install any python "package" as a command line tool:

This requires the `pyproject.toml` file included in this repo. Please see notes therein. 

From the terminal:

`cd scampiman`

`pip install .`

* This should install `scampiman` as a runnable command from the terminal

# Running `scampiman`

With a directory of unaligned `.bam` files from an ONT run:

```bash
scampiman -r proj1/bam_pass/barcode24 -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s barcode24 -o proj1_scampi -f bam -t directory --seqtech ont --temp $TMPDIR
```

With some unaligned `.bam` files from an ONT run:

```bash
scampiman -r proj1/bam_pass/barcode24/*bam -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s barcode24 -o proj1_scampi -f bam -t files --seqtech ont --temp $TMPDIR
```

With some `.fastq` files from an Illumina run:

```bash
scampiman -r my_fastqs/seq1.R1.fastq my_fastqs/seq1.R2.fastq -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s seq1 -o proj2_scampi -f fastq -t files --seqtech illumina --temp $TMPDIR
```


## Plotting data (not thoroughly tested/robust)

See conda environment requirements below.

This needs an index file in `.xlsx` format with (at least) the following header columns:
 - Barcode ID
 - Sample ID

```bash
Rscript scampiman/plot_script/plot_scampiman_batch1.R scampi_projects my_amplicons_projs1to4.pdf
```

# Other

`Dorado`, from ONT, is required but is only distributed as a binary. I used `v0.9.1`. Assumes `Dorado` in PATH.
The install instructions can be found [here](https://github.com/nanoporetech/dorado?tab=readme-ov-file#installation).

terminal commmand used to create `scampiman` conda environment:

```
$ conda create -n scampiman -c conda-forge -c bioconda python bioconda::minimap2 bioconda::samtools bioconda::pysam conda-forge::pandas 
```

terminal command to add `R` plotting capabilities 
```
$ conda activate scampiman

$ conda install -c conda-forge conda-forge::r-rprojroot conda-forge::r-tidyverse conda-forge::r-cowplot conda-forge::r-readxl
```