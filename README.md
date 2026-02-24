🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\
🌊🌊🌊🌊🌊🌊🌊🌊🦐🦐🦐🌊🦐🦐🌊🦐🦐🌊🦐🦐🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\
🌊🌊🌊🌊🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\
🌊🌊🦐🦐🌊🦐🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🦐🌊🦐🦐🌊🦐🦐🦐🌊🌊🌊🌊🌊🌊\
🌊🦐🌊🌊🌊🦐🌊🌊🦐🦐🦐🌊🦐🌊🌊🌊🦐🌊🦐🦐🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🌊🦐🌊\
🌊🌊🦐🌊🌊🦐🌊🌊🦐🌊🦐🌊🦐🌊🌊🌊🦐🌊🦐🌊🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🦐🦐🌊🦐🦐🌊🦐🌊\
🌊🌊🌊🦐🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🌊🌊🦐🌊🦐🌊🌊🌊🌊🌊🦐🌊🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🦐🌊\
🌊🦐🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🦐🌊🦐🌊🌊🦐🌊\
🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\
🌊🌊🌊🌊🌊🌊🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊

# scampiman
A pipeline to align, quality control, and summarize tiled amplicon coverage (of a virus, probably) from sequencing reads.

**Rationale:** We noticed 1) that tiled amplicon data can come in many forms from many technologies, and 2) errors introduced in library prep can lead to sequencing artifacts that, if not handled properly, can cause issues with downstream analysis.

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
2) `pysam`: sort, ampliconclip, ampliconstats
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
Scampiman assumes that paired‑end reads were generated on an Illumina or similar platform.

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

Note: consider making an isolated environment (conda or venv) for `scampiman`.

**Easiest Way**

`, Simply `pip install` with the release tarball link.

`pip install https://github.com/tiszalab/scampiman/archive/refs/tags/v0.1.0.tar.gz`

**Alternative Methods**:

1. clone this repo or download and unpack release.

2. `pip` install scampiman

From the terminal:

`cd scampiman`

`pip install .`

* Either method should install `scampiman` as a runnable command from the terminal.

# Running `scampiman`

**Highly Recommended**: quality filter reads before running scampiman with e.g. `fastp` (short reads) or `fastplong` (long reads)!

### With a directory of unaligned `.bam` files from an ONT run:

```bash
scampiman -r proj1/bam_pass/barcode24 -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s barcode24 -o proj1_scampi -f bam -t directory -c single-end --seqtech ont
```

You can also specify multiple directories:
```bash
scampiman -r flowcell1/bam_pass/barcode24 flowcell2/bam_pass/barcode24 -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s barcode24 -o proj1_scampi -f bam -t directory -c single-end --seqtech ont
```

### With some unaligned `.bam` files from an ONT run:

```bash
scampiman -r proj1/bam_pass/barcode24/*bam -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s barcode24 -o proj1_scampi -f bam -t files -c single-end --seqtech ont
```

It's better to use `-t directory` if you are using all files in a directory.

### With `.fastq` files:

from a paired-end Illumina run:
```bash
scampiman -r my_fastqs/seq1.R1.fastq my_fastqs/seq1.R2.fastq -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s seq1 -o proj2_scampi -f fastq -t files -c paired-end --seqtech illumina
```

Single-end  (e.g. ONT) works too:
```bash
scampiman -r my_fastqs/seq1.ONT.fastq -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s seq1 -o proj2_scampi -f fastq -t files -c single-end --seqtech ont
```

### Keeping aligned `.bam` files for downstream analysis:

Commonly, you will want to keep the properly filtered `.bam` file to run downstream analysis to determine lineage, derive consensus genome, or analyze allele frequency.

```bash
scampiman -r my_fastqs/seq1.ONT.fastq -b SARS-CoV-2.ARTIC_5.3.2.primer.bed -g sars_cov2_MN908947.3.fasta -s seq1 -o proj2_scampi -f fastq -t files -c single-end --seqtech ont --keep bam
```


## Plotting data (not thoroughly tested/robust)

See conda environment requirements below.

This needs an index file in `.xlsx` format with (at least) the following header columns:
 - Barcode ID
 - Sample ID

```bash
Rscript scampiman/plot_script/plot_scampiman_batch1.R scampi_projects my_amplicons_projs1to4.pdf
```


terminal command to add `R` plotting capabilities 
```
conda activate scampiman

conda install -c conda-forge conda-forge::r-rprojroot conda-forge::r-tidyverse conda-forge::r-cowplot conda-forge::r-readxl
```


