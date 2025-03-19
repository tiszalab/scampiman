# scampiman
pipeline to align and summarize tiled amplicon coverage from raw sequences

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