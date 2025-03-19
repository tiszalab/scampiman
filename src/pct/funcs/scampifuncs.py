import logging
import argparse
from subprocess import Popen

# since I'm logging, I need the same logger as in the main script
logger = logging.getLogger("pct_logger")

def cat_bams_dir(readdir: str, cpus: int, outf: str) -> Popen:

    return subprocess.Popen(['samtools', 'cat', '--threads', cpus,
                    '-o', outf,
                    f'{readdir}/*bam'],
                    stdout=PIPE, stderr=STDOUT)

def cat_bams_files(readfiles: str, cpus: int, outf: str) -> Popen:

    return subprocess.Popen(['samtools', 'cat', '--threads', cpus,
                    '-o', outf,
                    readfiles],
                    stdout=PIPE, stderr=STDOUT)

def dorado_al(bam: str, cpus: int, outf: str, ref: str):

    sortbam = open(outf, 'a')

    dorado_command = ['dorado', 'aligner', '-t', cpus, ref, bam]

    # Second command-line
    view_command = ['samtools', 'view', '-F', '4', '-Sb']
    
    # Third command-line
    sort_command = ['samtools', 'sort', '-o', outf]

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

    return output.decode()

def mini2_al(fq: str, cpus: int, outf: str, ref: str):

    sortbam = open(outf, 'a')

    mini_command = ['minimap2', '-ax', 'sr', '-t', cpus, ref, fq]

    # Second command-line
    view_command = ['samtools', 'view', '-F', '4', '-Sb']
    
    # Third command-line
    sort_command = ['samtools', 'sort', '-o', outf]

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

    return output.decode()

def ampclip(bed: str, sortbam: str, outf: str):

    return subprocess.Popen(['samtools', 'ampliconclip', '-b', bed,
                    '-o', outf,
                    sortbam],
                    stdout=PIPE, stderr=STDOUT)

def ampstats(bed: str, clipbam: str, outf: str):
    
    return subprocess.Popen(['samtools', 'ampliconstats', 
                    '-o', outf,
                    bed, clipbam
                    ],
                    stdout=PIPE, stderr=STDOUT)


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

