#!/usr/bin/env python

import argparse
import sys, os
import pysam
import shutil
import time
import logging
from datetime import timedelta

# importing this packages functions from funcs/scampifuncs.py
# the try/except loop allows this tool to be run as `python scampiman.py` (.funcs), for development
# or be run as a true command line tool, `scampiman`, after being installed with pip (funcs)
try:
    from .funcs import scampifuncs as scaf
except:
    from funcs import scampifuncs as scaf


def scampiman():

    starttime = time.perf_counter()

    # I like to keep this print statement for debugging
    print(f"this script dir: {os.path.dirname(__file__) }")

    __version__ = "0.1.3"

    toolname = "scampiman"

    # this variable isn't used in this template, but it's useful for almost all bioinformatics
    def_workdir = os.getcwd()

    parser = argparse.ArgumentParser(prog=f'{toolname}',
                                     description=f'{toolname}. amplicon tool. Visit '
                                    f'https://github.com/tiszalab/{toolname} for help.\n'
                                    f'Version {str(__version__)}')
    
    parser.add_argument("-r", "--reads", nargs="+",
                        dest="READS", required=True, 
                        help='read file(s) or directory containing read files. \
                        fastq or unaligned bam format'
                        )
    parser.add_argument("-b", "--bed", type=str,
                        dest="bed", required=True, 
                        help='bed file of primers'
                        )
    parser.add_argument("-g", "--genome", type=str,
                        dest="genome", required=True, 
                        help='genome file. Accession/sequence must match bed file'
                        )
    parser.add_argument("-s", "--sample", 
                        dest="SAMPLE", type=str, required=True, 
                        help='Sample name. No space characters, please.'
                        )
    parser.add_argument("-o", "--output_dir", 
                        dest="OUTPUT_DIR", type=str, required=True, 
                        help='Output directory name. Will be created if it does not exist. \
                        Can be shared with other samples. No space characters, please. '
                        )
    parser.add_argument("-f", "--read-format",
                        dest="rfmt", default="fastq", 
                        choices=['fastq', 'bam'],
                        help='read file format'
                        )
    parser.add_argument("-c", "--read-configuration",
                        dest="rcon", default="single-end", 
                        choices=['single-end', 'paired-end'],
                        help='read file configuration'
                        )
    parser.add_argument("-t", "--input-type",
                        dest="intype", default="files", 
                        choices=['files', 'directory'],
                        help='read file format'
                        )
    parser.add_argument("-@", "--cpus",
                        dest="cpus", type=int,
                        help='number of cpus to use'
                        )
    parser.add_argument("--seqtech", dest="SEQTECH",
                        default='ont', 
                        choices=['illumina', 'ont'],
                        help='Which sequencing technology produced the reads?'
                        )
    parser.add_argument("--temp", 
                        dest="TEMP_DIR", type=str, default='default',
                        help='path of temporary directory. Default is {OUTPUT_DIR}/{SAMPLE}_temp/')
    parser.add_argument("--keep", 
                        dest="KEEP", type=scaf.str2bool, default='False',
                        help='True of False. Keep the intermediate files, located in the temporary directory? These can add up, so it is not recommended if space is a concern.') 
    parser.add_argument("-wd", "--working_directory", dest="c_workdir", type=str, default=def_workdir, 
                        help=f"Default: {def_workdir} -- \
                        Set working directory with absolute or relative path. \
                        Run directory will be created within.")
    
    args = parser.parse_args()


    #### define logger #####
    logger = logging.getLogger("pct_logger")
    logger.setLevel(logging.DEBUG)
    # stream gets printed to terminal
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.ERROR)

    # file gets saved to a specified file
    file_handler = logging.FileHandler(os.path.join(def_workdir, 
                                                    f"scampi{args.SAMPLE}.log"))
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    #########################
    # get console width    
    scaf.shrimp_header(__version__)
    
    ## make directories
    out_directory = os.path.join(str(args.c_workdir), str(args.OUTPUT_DIR))

    if not args.cpus:
        def_CPUs = os.cpu_count()
    else:
        def_CPUs = args.cpus

    if not os.path.isdir(out_directory):
        os.makedirs(out_directory)

    if str(args.TEMP_DIR) == 'default':
        sca_temp = os.path.join(out_directory, f'{str(args.SAMPLE)}_temp')
    else:
        sca_temp = str(args.TEMP_DIR)

    if str(sca_temp) == str(args.OUTPUT_DIR):
        logger.error("temporary directory and output directory are the same.")
        logger.error("This is not allowed. Exiting.")
        sys.exit()

    if not os.path.isdir(sca_temp):
        os.makedirs(sca_temp)
        logger.info(f"temp dir: {sca_temp}")
        logger.info(os.listdir(sca_temp))


    #check if files exist
    for inf in [args.genome, args.bed]:
        if not os.path.isfile(inf):
            logger.error(f'input file not found at {inf}. exiting.')
            sys.exit()
    
    #check if bed and genome files have same accessions
    if not scaf.bed_genome_match(
        str(args.bed),
        str(args.genome)
    ):
        logger.error("accessions in bed file and genome file do not match.")
        logger.error("Exiting.")
        sys.exit()
    
    # if args.rfmt is "bam" and args.rcon is "paired-end", exit
    if str(args.rfmt) == "bam" and str(args.rcon) == "paired-end":
        logger.error("bam format is not supported for paired-end reads.")
        logger.error("Exiting.")
        sys.exit()

    #READ_STR = ' '.join(map(str,args.READS))
    logger.info(f'read string: ')
    logger.info(f'{args.READS}')

    align_starttime = time.perf_counter()
    logger.info(f'Starting alignment.')
    
    # run main pipeline
    try:
        if str(args.rfmt) == "bam":
            logger.info(f'> bam option ')
            print(f"🦐 bam option ")
        elif str(args.rfmt) == "fastq":
            logger.info(f'> fastq option ')
            print(f"🦐 fastq option ")
        if str(args.intype) == "files":
            logger.info(f'> files option ')
            print(f"🦐 files option ")
        elif str(args.intype) == "directory":
            logger.info(f'> directory option ')
            print(f"🦐 directory option ")

        if str(args.rcon) == "single-end":
            logger.info(f'> single-end option ')
            print(f"🦐 single-end reads ")
            print(f"✨✨ Let's go! ✨✨")
            try:
                scaf.shrimp_progress(1, 0, 0, "preprocessing")
                reads_list = scaf.file_paths(args.READS, args.rfmt, args.rcon, args.intype)

                alignstats = scaf.mappy_al_single(
                    args.rfmt,
                    def_CPUs,
                    args.SEQTECH,
                    str(args.genome),
                    os.path.join(sca_temp, f'{str(args.SAMPLE)}.bam'),
                    os.path.join(sca_temp, f'{str(args.SAMPLE)}.failed.bam'),
                    reads_list
                )

                scaf.stats_tsv(alignstats, 1, os.path.join(args.OUTPUT_DIR, f'{str(args.SAMPLE)}.summarystats.tsv'))
            except:
                logger.error("Failed to align single-end read files.")

        elif str(args.rcon) == "paired-end":
            logger.info(f'> paired-end option ')
            print(f"🦐 paired-end reads ")
            print(f"✨✨ Let's go ✨✨")

            try:
                scaf.shrimp_progress(2, 0, 0, "preprocessing")
                read1_list, read2_list = scaf.file_paths(args.READS, args.rfmt, args.rcon, args.intype)

                alignstats = scaf.mappy_al_paired(
                    def_CPUs,
                    args.SEQTECH,
                    str(args.genome),
                    os.path.join(sca_temp, f'{str(args.SAMPLE)}.bam'),
                    os.path.join(sca_temp, f'{str(args.SAMPLE)}.failed.bam'),
                    read1_list,
                    read2_list
                )

                scaf.stats_tsv(alignstats, 1, os.path.join(args.OUTPUT_DIR, f'{str(args.SAMPLE)}.summarystats.tsv'))
            except:
                logger.error("Failed to align paired-end read files.")

    except Exception as e:
        logger.error("Read Alignment ERROR: ")
        logger.error(e)

    align_endtime = time.perf_counter()

    time_taken = align_endtime - align_starttime

    logger.info(f"> alignment step took {timedelta(seconds=time_taken)}")

    logger.info(f'Starting amplicon analysis.')
    
    amp_starttime = time.perf_counter()

    if os.path.isfile(os.path.join(sca_temp, f'{str(args.SAMPLE)}.sort.bam')):
        try:
            scaf.shrimp_progress(3, 0, 0, "amp")

            logger.info(f'samcov')
            pysam.samtools.coverage(
                '-o', os.path.join(
                    args.OUTPUT_DIR, 
                    f'{str(args.SAMPLE)}.samcov.tsv'
                    ),
                os.path.join(sca_temp, f'{str(args.SAMPLE)}.sort.bam') 
            )
            amp_endtime = time.perf_counter()
            time_taken = amp_endtime - amp_starttime
            scaf.shrimp_progress(3, 1, time_taken, "amp")
            logger.info(f'ampclip')
            logger.info(f'longest amplicon length: {int(scaf.ampliconstats_length(args.bed)/2)}bp')

            pysam.samtools.ampliconclip(
                '--both-ends',
                '-@', str(def_CPUs),
                '-b', args.bed,
                '-o', 
                os.path.join(
                    sca_temp, 
                    f'{str(args.SAMPLE)}.ampclip.bam'
                    ),
                os.path.join(sca_temp, f'{str(args.SAMPLE)}.sort.bam')
            )
            amp_endtime = time.perf_counter()
            time_taken = amp_endtime - amp_starttime
            scaf.shrimp_progress(3, 2, time_taken, "amp")
            pysam.samtools.ampliconstats(
                '-@', str(def_CPUs),
                '-l', str(scaf.ampliconstats_length(args.bed)),
                '-o', 
                os.path.join(
                    args.OUTPUT_DIR,
                    f'{str(args.SAMPLE)}.ampliconstats.tsv'
                ),
                args.bed,
                os.path.join(
                    sca_temp, 
                    f'{str(args.SAMPLE)}.ampclip.bam'
                    )
            )

            logger.info(f"temp dir: {sca_temp}")
            logger.info(os.listdir(sca_temp))


            logger.info(f'amptable')
            scaf.amptable(
                os.path.join(
                    args.OUTPUT_DIR,
                    f'{str(args.SAMPLE)}.ampliconstats.tsv'
                ),
                args.SAMPLE,
                os.path.join(
                    args.OUTPUT_DIR,
                    f'{str(args.SAMPLE)}.amplicontable.tsv'
                )
            )

            amp_endtime = time.perf_counter()
            time_taken = amp_endtime - amp_starttime            
            scaf.shrimp_progress(3, 3, time_taken, "amp")

        except Exception as e:
            logger.error("Amplicon Analysis ERROR: ")
            logger.error(e)

    logger.info(f"### scampiman outputs: ")
    for fin in [
        f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.summarystats.tsv',
        f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.amplicontable.tsv',
        f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.ampliconstats.tsv',
        f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.samcov.tsv'
    ]:
        if os.path.isfile(fin):
            logger.info(f"### detected - {fin}")
            print(f"  🍤  {fin}")
        else:
            logger.info(f"!!!! Not found - {fin}")
            print(f"  Ⓧ NOT FOUND - {fin}")
    
    if os.path.isdir(sca_temp) and not args.KEEP:
        logger.info(f"removing temp files in: {sca_temp}")
        shutil.rmtree(sca_temp)

    endtime = time.perf_counter()

    time_taken = endtime - starttime

    logger.info(f"> scampiman took {timedelta(seconds=time_taken)}")

# this is how to run the pct function by calling this script from the command line
if __name__ == "__main__":
    scampiman()

