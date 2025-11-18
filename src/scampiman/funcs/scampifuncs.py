import logging
import argparse
import pysam
import sys, os
import pandas as pd
import subprocess
from subprocess import Popen, PIPE, STDOUT

# since I'm logging, I need the same logger as in the main script
logger = logging.getLogger("pct_logger")

def pysam_cat(reads: str, ftype: str, outf: str):

    if ftype == "directory":
        bam_list = []
        for bam in os.listdir(reads):
            if bam.endswith('.bam'):
                f = os.path.join(reads, bam)

                if os.path.isfile(f) and os.path.getsize(f) > 0:
                    bam_list.append(f)

    elif ftype == "files":
        bam_list = []
        for bam in reads.split():
            if bam.endswith('.bam'):
                f = os.path.abspath(bam)

                if os.path.isfile(f) and os.path.getsize(f) > 0:
                    bam_list.append(f)

    pysam.samtools.cat(
        '-o', outf,
        *bam_list
    )

def dorado_al(bam: str, cpus: int, outf: str, ref: str):

    sortbam = open(outf, 'a')

    dorado_command = ['dorado', 'aligner', '-t', str(cpus), ref, bam]

    # Second command-line
    view_command = ['samtools', 'view', '-Sb']
    logger.info(view_command)
    
    # Third command-line
    sort_command = ['samtools', 'sort', '-@', str(cpus), '-o', outf]

    # Launch first process
    dorado_process = subprocess.Popen(
        dorado_command,
        stdout=subprocess.PIPE
        )

    # Launch second process and connect it to the first one
    view_process = subprocess.Popen(
        view_command, 
        stdin=dorado_process.stdout, 
        stdout=subprocess.PIPE
    )

    # Launch third process and connect it to the first one
    sort_process = subprocess.Popen(
        sort_command, 
        stdin=view_process.stdout, 
        stdout=sortbam
    )

    # Let stream flow between them
    output, _ = sort_process.communicate()
    if output:
        return output
    else:
        return None


def mini2_al(fq: list, cpus: int, outf: str, ref: str, tech: str):
    print("mini2_al")
    print(fq)
    sortbam = open(outf, 'a')

    if tech == "ont":
        targ = "lr:hq"
    elif tech == "illumina":
        targ = "sr"

    fq_list = []
    for faq in fq:
        if os.path.isfile(faq) and os.path.getsize(faq) > 0:
            fq_list.append(os.path.abspath(faq))
    mini_command = ['minimap2', '-ax', targ, '-t', str(cpus), ref, *fq_list]

    # Second command-line
    view_command = ['samtools', 'view', '-Sb']
    
    # Third command-line
    sort_command = ['samtools', 'sort', '-@', str(cpus), '-o', outf]

    # Launch first process
    mini_process = subprocess.Popen(
        mini_command,
        stdout=subprocess.PIPE
        )

    # Launch second process and connect it to the first one
    view_process = subprocess.Popen(
        view_command, 
        stdin=mini_process.stdout, 
        stdout=subprocess.PIPE
    )

    # Launch third process and connect it to the first one
    sort_process = subprocess.Popen(
        sort_command, 
        stdin=view_process.stdout, 
        stdout=sortbam
    )

    # Let stream flow between them
    output, _ = sort_process.communicate()

    if output:
        return output
    else:
        return None


def amptable(ampstats: str, sampid: str) -> pd.DataFrame:

    odf = pd.DataFrame(columns=[
        'accession',
        'amplicon_number',
        'lprimer',
        'rprimer'
        ]
    )
    with open(ampstats, 'r') as sam:
        for line in sam:
            line = line.strip()
            if line.startswith("AMPLICON"):
                ampinfo = line.split('\t')[1:]
                odf.loc[len(odf)] = ampinfo
            if line.startswith("FREADS"):
                columns = line.split('\t')[2:]
                if len(odf) == len(columns):
                    odf['amplicon_reads'] = columns
                else:
                    raise Exception("amplicon_reads data incorrect length")

            if line.startswith("FVDEPTH"):
                columns = line.split('\t')[2:]
                if len(odf) == len(columns):
                    odf['full_length_depth'] = columns
                else:
                    raise Exception("full_length_depth data incorrect length")
            if line.startswith("FDEPTH"):
                columns = line.split('\t')[2:]
                if len(odf) == len(columns):
                    odf['avg_depth'] = columns
                else:
                    raise Exception("avg_depth data incorrect length")
            if line.startswith("FPCOV-1"):
                columns = line.split('\t')[2:]
                if len(odf) == len(columns):
                    odf['amplicon_coverage'] = columns
                else:
                    raise Exception("amplicon_coverage data incorrect length")
    odf['sample_ID'] = sampid
    return odf

def bed_genome_match(bed: str, genome: str) -> bool:
    bedlist = []
    with open(bed, 'r') as bf:
        for line in bf:
            line = line.strip()
            acc = line.split('\t')[0]
            if acc not in bedlist:
                bedlist.append(acc)
    
    genomelist = []
    with open(genome, 'r') as gf:
        for line in gf:
            line = line.strip()
            if line.startswith(">"):
                acc = line.lstrip(">").split(' ')[0]
                if acc not in genomelist:
                    genomelist.append(acc)

    if set(bedlist) == set(genomelist):
        return True
    else:
        logger.info(f'BED: {set(bedlist)}')
        logger.info(f'GENOME: {set(genomelist)}')
        return False


def str2bool(v):
    '''
    converts bool-like strings to True/False boolean objects
    '''
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def ampliconstats_length(bed: str) -> int:
    amplicons = {}
    with open(bed, 'r') as bf:
        for line in bf:
            stripped = line.strip()
            columns = stripped.split('\t')
            primer_name = columns[3]
            if '_LEFT' in primer_name.upper():
                amplicon_id = primer_name.rsplit('_LEFT', maxsplit=1)[0]
                left_start = int(columns[1])
                amplicons.setdefault(amplicon_id, {})['left_start'] = left_start
            elif '_RIGHT' in primer_name.upper():
                amplicon_id = primer_name.rsplit('_RIGHT', maxsplit=1)[0]
                right_end = int(columns[2])
                amplicons.setdefault(amplicon_id, {})['right_end'] = right_end
            else:
                raise ValueError("primer not left or right")
    longest = None
    max_length = -1
    for amp_id, coords in amplicons.items():
        if 'left_start' not in coords or 'right_end' not in coords:
            continue
        length = coords['right_end'] - coords['left_start']
        if length > max_length:
            max_length = length
            longest = {
                'amplicon_number': amp_id,
                'left_start': coords['left_start'],
                'right_end': coords['right_end'],
                'length': length
            }
    if longest is None:
        raise ValueError('No complete amplicon pairs found in BED file')
    max_length += 50
    return max_length