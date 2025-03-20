#!/usr/bin/env python

import argparse
import sys, os
#import subprocess
import shutil
import subprocess
from subprocess import Popen, PIPE, STDOUT
from pathlib import Path
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

    __version__ = "0.0.1"

    toolname = "scampiman"

    # this variable isn't used in this template, but it's useful for almost all bioinformatics
    def_CPUs = os.cpu_count()

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
    parser.add_argument("-t", "--input-type",
                        dest="intype", default="files", 
                        choices=['files', 'directory'],
                        help='read file format'
                        )

    parser.add_argument("--seqtech", dest="SEQTECH",
                        default='illumina', 
                        choices=['illumina', 'ont'],
                        help='Which sequencing technology produced the reads?'
                        )
    parser.add_argument("--temp", 
                        dest="TEMP_DIR", type=str, default='default',
                        help='path of temporary directory. Default is {OUTPUT_DIR}/{SAMPLE}_temp/')   
    parser.add_argument("-wd", "--working_directory", dest="c_workdir", type=str, default=def_workdir, 
                        help=f"Default: {def_workdir} -- \
                        Set working directory with absolute or relative path. \
                        Run directory will be created within.")
    # you need to run this line to parse the arguments
    args = parser.parse_args()


    #### define logger #####
    logger = logging.getLogger("pct_logger")
    logger.setLevel(logging.DEBUG)
    # stream gets printed to terminal
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)

    # file gets saved to a specified file
    file_handler = logging.FileHandler(os.path.join(def_workdir, 
                                                    f"scampi{args.SAMPLE}.log"))
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    #########################

    ## make directories
    out_directory = os.path.join(str(args.c_workdir), str(args.OUTPUT_DIR))

    if not os.path.isdir(out_directory):
        os.makedirs(out_directory)

    samp_out_dir = os.path.join(out_directory, str(args.SAMPLE))

    if not os.path.isdir(samp_out_dir):
        os.makedirs(samp_out_dir)

    if str(args.TEMP_DIR) == 'default':
        sca_temp = os.path.join(samp_out_dir, f'{str(args.SAMPLE)}_temp')
    else:
        sca_temp = str(args.TEMP_DIR)

    if str(sca_temp) == str(args.OUTPUT_DIR):
        logger.error("temporary directory and output directory are the same.")
        logger.error("This is not allowed. Exiting.")
        sys.exit()

    if not os.path.isdir(sca_temp):
        os.makedirs(sca_temp)


    #check if files exist
    for inf in [args.genome, args.bed]:
        if not os.path.isfile(inf):
            logger.error(f'input file not found at {inf}. exiting.')
            sys.exit()

    align_starttime = time.perf_counter()

    logger.info(f'Starting alignment.')

    READ_STR = ' '.join(map(str,args.READS))
    logger.info(f'read string: ')
    logger.info(f'{READ_STR}')

    try:
        if str(args.rfmt) == "bam":
            logger.info(f'> bam option ')
            if str(args.intype) == "files":
                logger.info(f'> files option ')

                if len(args.READS) > 1:
                    scaf.cat_bams_files(
                        READ_STR, 
                        def_CPUs, 
                        f'{sca_temp}/{str(args.SAMPLE)}_cat.bam'
                    )
                    scaf.dorado_al(
                        f'{sca_temp}/{str(args.SAMPLE)}_cat.bam', 
                        def_CPUs, 
                        f'{sca_temp}/{str(args.SAMPLE)}.sort.bam', 
                        str(args.genome)
                    )
                else:
                    scaf.dorado_al(
                        READ_STR, 
                        def_CPUs, 
                        f'{sca_temp}/{str(args.SAMPLE)}.sort.bam', 
                        str(args.genome)
                    )
            if str(args.intype) == "directory":
                logger.info(f'> directory option ')
                try: 
                    scaf.cat_bams_dir(
                        READ_STR,
                        def_CPUs,
                        f'{sca_temp}/{str(args.SAMPLE)}_cat.bam'
                    )
                except Exception as e:
                    logger.error("samtools cat not successful:")
                    logger.error(e)


                scaf.dorado_al(
                    f'{sca_temp}/{str(args.SAMPLE)}_cat.bam', 
                    def_CPUs, 
                    f'{sca_temp}/{str(args.SAMPLE)}.sort.bam', 
                    str(args.genome)
                )
        if str(args.rfmt) == "fastq":
            logger.info(f'> fastq option ')
            scaf.mini2_al(
                READ_STR,
                def_CPUs,
                f'{sca_temp}/{str(args.SAMPLE)}.sort.bam', 
                str(args.genome)
            )
    except Exception as e:
        logger.error("Read Alignment ERROR: ")
        logger.error(e)

    align_endtime = time.perf_counter()

    time_taken = align_endtime - align_starttime

    logger.info(f"> alignment step took {timedelta(seconds=time_taken)}")

    logger.info(f'Starting amplicon analysis.')

    logger.info(f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.ampliconstats.tsv')
    logger.info(str(args.bed))
    logger.info(f'{sca_temp}/{str(args.SAMPLE)}.ampclip.bam')

    if os.path.isfile(f'{sca_temp}/{str(args.SAMPLE)}.sort.bam'):
        try:
            logger.info(f'ampclip')
            scaf.ampclip(
                str(args.bed),
                f'{sca_temp}/{str(args.SAMPLE)}.sort.bam',
                f'{sca_temp}/{str(args.SAMPLE)}.ampclip.bam'
            )
            ##### SOMETHING IS GOING ON WITH PYTHON RECOGNIZING f'{sca_temp}/{str(args.SAMPLE)}.ampclip.bam' ####
            if os.path.exists(f'{sca_temp}/{str(args.SAMPLE)}.sort.bam'):
                logger.info(f'{sca_temp}/{str(args.SAMPLE)}.sort.bam is found')
            else:
                logger.info(f'{sca_temp}/{str(args.SAMPLE)}.sort.bam NOT found')

            if os.path.exists(f'{sca_temp}/{str(args.SAMPLE)}.ampclip.bam'):
                logger.info(f'{sca_temp}/{str(args.SAMPLE)}.ampclip.bam is found')
            else:
                logger.info(f'{sca_temp}/{str(args.SAMPLE)}.ampclip.bam NOT found')
                
            logger.info(f'ampstats')
            scaf.ampstats(
                str(args.bed),
                f'{sca_temp}/{str(args.SAMPLE)}.sort.bam',
                f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.ampliconstats.tsv'
            )
#            clipf = os.path.join(sca_temp, args.SAMPLE + '.ampclip.bam')
#            scaf.run_samtools_ampliconstats(
#                str(args.bed), 
#                #f'{sca_temp}/{str(args.SAMPLE)}.ampclip.bam', 
#                os.path.join(sca_temp, args.SAMPLE + '.sort.bam'),
#                #str(clipf),
#                f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.ampliconstats.tsv'
#                
#            )
#            ampf = os.path.join({str(args.OUTPUT_DIR)}, f'{str(args.SAMPLE)}.ampliconstats.tsv')
#            logger.info(f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.ampliconstats.tsv')
#            try:
#                subprocess.Popen(['samtools', 'ampliconstats', 
#                        '-o', 
#                        f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.ampliconstats.tsv',
#                        str(args.bed), 
#                        f'{sca_temp}/{str(args.SAMPLE)}.ampclip.bam'
#                        ], 
#                        shell=False
#                )
#            except Exception as e:
#                logger.error(e)

            logger.info(f'samcov')
            scaf.samcov(
                f'{sca_temp}/{str(args.SAMPLE)}.sort.bam',
                f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.samcov.tsv'
            )
            logger.info(f'amptable')
            odf = scaf.amptable(
                f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.ampliconstats.tsv'
            )

            odf.to_csv(
                f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.amplicontable.tsv',
                sep = "\t",
                index=False
            )

        except Exception as e:
            logger.error("Amplicon Analysis ERROR: ")
            logger.error(e)

    logger.info(f"### scampiman outputs: ")
    for fin in [
        f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.amplicontable.tsv',
        f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.ampliconstats.tsv',
        f'{str(args.OUTPUT_DIR)}/{str(args.SAMPLE)}.samcov.tsv'
    ]:
        if os.path.isfile(fin):
            logger.info(f"### detected - {fin}")
        else:
            logger.info(f"!!!! Not found - {fin}")
    
#    if os.path.isdir(sca_temp):
#        shutil.rmtree(sca_temp)

    endtime = time.perf_counter()

    time_taken = endtime - starttime

    logger.info(f"> scampiman took {timedelta(seconds=time_taken)}")

# this is how to run the pct function by calling this script from the command line
if __name__ == "__main__":
    scampiman()

