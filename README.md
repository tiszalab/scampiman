🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊
🌊🌊🌊🌊🌊🌊🌊🌊🦐🦐🦐🌊🦐🦐🌊🦐🦐🌊🦐🦐🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊
🌊🌊🌊🌊🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊
🌊🌊🦐🦐🌊🦐🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🦐🌊🦐🦐🌊🦐🦐🦐🌊🌊🌊🌊🌊🌊
🌊🦐🌊🌊🌊🦐🌊🌊🦐🦐🦐🌊🦐🌊🌊🌊🦐🌊🦐🦐🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🌊🦐🌊
🌊🌊🦐🌊🌊🦐🌊🌊🦐🌊🦐🌊🦐🌊🌊🌊🦐🌊🦐🌊🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🦐🦐🌊🦐🦐🌊🦐🌊
🌊🌊🌊🦐🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🌊🌊🦐🌊🦐🌊🌊🌊🌊🌊🦐🌊🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🦐🌊
🌊🦐🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🦐🌊🦐🌊🌊🦐🌊
🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊
🌊🌊🌊🌊🌊🌊🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊

# scampiman
A pipeline to align and summarize tiled amplicon coverage (of a virus, probably) from raw sequences.

**Requires: Reads, Reference Genome(s), Primer `.bed` File**

**Produces: Alignment Summary, Samtools Ampliconstats File, Table of Amplicon Coverage `.tsv`**

 - Input formats
    `.fastq` or `.bam`
 - Input logic
    `file(s)` or `directory` (with files)
 - Seq tech
    `illumina short read` or `ONT`
 - Read config
    `single-end` or `paired-end`

1) Align reads to reference (`mappy`) and filter unwanted alignments
2) pysam: sort, ampliconclip, ampliconstats
3) Parse ampliconstats output into table, output `.tsv`

# Alignment Filtering

Subpar alignments are filtered out before amplicon analysis is performed. This step attempts to remove issues that may have arisen during library preparation, for either single‑ or paired‑end reads, that can cause misrepresentation of amplicon diversity.

The number of removed alignments is reported in the alignment summary as 'removed_reads_primary' and is saved within the failed.bam output along with unmapped reads.

## Single‑end Reads
The filtering parameters for single‑end reads are designed to correct for ligation‑based errors that may occur, particularly in ONT ligation‑based sequencing kits (e.g., SQK‑NBD114).

The filtering parameters are as follows:
- Removes reads with supplementary alignments that overlap <50% with the primary alignment’s reference region.
- Removes reads that produce supplementary alignments mapping to the same strand as the primary alignment.

Removal of these reads is important because it accounts for:
1. ligation between amplicons originating from different regions of the genome.
2. ligation between segments originating from different sources. For example, different barcodes of ONT kits.

## Paired‑end Reads
Scampiman assumes that paired‑end reads were generated on an Illumina platform.

The filtering parameters are as follows:
- Removes paired reads that align to the same strand.
- Removes paired reads that have unequal numbers of alignments (indicating mapping error).
- Removes paired reads where one mate is unmapped.
- Removes paired reads whose reference alignments do not overlap.

Removal of these reads is important because it accounts for:
1. Illumina's platform sequencing  paired reads from opposing strands of the same DNA fragment.
2. potential ligation or mapping errors.
3. the necessity for the entire (gap-less) amplicon to be represented in the analysis.

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
scampiman -r proj1/bam_pass/barcode24 -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s barcode24 -o proj1_scampi -f bam -t directory -c single-end --seqtech ont  --temp $TMPDIR
```

With some unaligned `.bam` files from an ONT run:

```bash
scampiman -r proj1/bam_pass/barcode24/*bam -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s barcode24 -o proj1_scampi -f bam -t files -c single-end --seqtech ont --temp $TMPDIR
```

With some `.fastq` files from an Illumina run:

```bash
scampiman -r my_fastqs/seq1.R1.fastq my_fastqs/seq1.R2.fastq -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s seq1 -o proj2_scampi -f fastq -t files -c paired-end --seqtech illumina --temp $TMPDIR
```


## Plotting data (not thoroughly tested/robust)

See conda environment requirements below.

This needs an index file in `.xlsx` format with (at least) the following header columns:
 - Barcode ID
 - Sample ID

```bash
Rscript scampiman/plot_script/plot_scampiman_batch1.R scampi_projects my_amplicons_projs1to4.pdf
```


terminal commmand used to create `scampiman` conda environment:

```
$ conda create -n scampiman -c conda-forge -c bioconda python bioconda::mappy bioconda::pysam conda-forge::pandas 
```

terminal command to add `R` plotting capabilities 
```
$ conda activate scampiman

$ conda install -c conda-forge conda-forge::r-rprojroot conda-forge::r-tidyverse conda-forge::r-cowplot conda-forge::r-readxl
```


