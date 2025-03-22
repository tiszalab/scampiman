import logging
import argparse
import pysam
import sys, os
import pandas as pd
import subprocess
from subprocess import Popen, PIPE, STDOUT

# since I'm logging, I need the same logger as in the main script
logger = logging.getLogger("pct_logger")

def cat_bams_dir(readdir: str, cpus: int, outf: str) -> Popen:
    bam_list = []
    for bam in os.listdir(readdir):
        if bam.endswith('.bam'):
            f = os.path.join(readdir, bam)

            if os.path.isfile(f) and os.path.getsize(f) > 0:
                bam_list.append(f)

    #logger.info(' '.join(bam_list))
    bamfs: str = ' '.join(bam_list)

    return subprocess.Popen(['samtools', 'cat',
                    '-o', outf,
                    bamfs
                    ],
                    stdout=PIPE, stderr=STDOUT)

def pysam_cat(reads: str, ftype: str, outf: str):

    if ftype == "directory":
        bam_list = []
        for bam in os.listdir(reads):
            if bam.endswith('.bam'):
                f = os.path.join(reads, bam)

                if os.path.isfile(f) and os.path.getsize(f) > 0:
                    bam_list.append(f)
        bamfs: str = ' '.join(bam_list)

        pysam.cat(
            bam_list,
            save_stdout = outf
        )
    elif ftype == "files":

        pysam.cat(
            reads,
            save_stdout = outf
        )  

def cat_bams_files(readfiles: str, cpus: int, outf: str) -> Popen:

    return subprocess.Popen(['samtools', 'cat',
                    '-o', outf,
                    readfiles],
                    stdout=PIPE, stderr=STDOUT)

def dorado_al(bam: str, cpus: int, outf: str, ref: str):

    sortbam = open(outf, 'a')


    dorado_command = ['dorado', 'aligner', '-t', str(cpus), ref, bam]

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
    if output:
        return output
    else:
        return None


def mini2_al(fq: str, cpus: int, outf: str, ref: str):

    sortbam = open(outf, 'a')

    mini_command = ['minimap2', '-ax', 'sr', '-t', str(cpus), ref, fq]

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

    if output:
        return output
    else:
        return None

def ampclip(bed: str, sortbam: str, outf: str):

    return subprocess.Popen(['samtools', 'ampliconclip', '-b', bed,
                    '-o', outf,
                    sortbam],
                    stdout=PIPE, stderr=STDOUT)

def ampstats(bed: str, clipbam: str, outf: str):

#    ampout = subprocess.run(['samtools', 'ampliconstats', 
#                    '-o', outf,
#                    bed, clipbam
#                    ],
#                    capture_output=True, text=True, check=True)
    
#    return ampout.stderr
#    ampf = open(outf, 'a')
    return subprocess.Popen(['samtools', 'ampliconstats', 
                    '-o', outf,
                    bed, clipbam
                    ],
                    stdout=PIPE, stderr=STDOUT)



def run_samtools_ampliconstats(bed_file, clip_bam_file, output_file):
    """
    Run samtools ampliconstats with the specified parameters.
    
    Args:
        bed_file (str): Path to the BED file with amplicon regions
        clip_bam_file (str): Path to the clipped BAM file
        output_file (str): Path to write the output results
        
    Returns:
        bool: True if the command executed successfully, False otherwise
    """
    try:
        # Ensure the input files exist
        if not os.path.isfile(bed_file):
            raise FileNotFoundError(f"BED file not found: {bed_file}")
        
        if not os.path.isfile(clip_bam_file):
            raise FileNotFoundError(f"BAM file not found: {clip_bam_file}")
        
        # Build the command
        command = ["samtools", "ampliconstats", "-o", output_file, bed_file, clip_bam_file]
        
        # Run the command
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        print(f"Command completed successfully. Output saved to: {output_file}")
        return True
        
    except FileNotFoundError as e:
        print(f"Error: {str(e)}")
        return False
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        return False


def samcov(sortbam: str, outf: str):

    covf = open(outf, 'a')
    return subprocess.Popen(['samtools', 'coverage', sortbam],
                            stdout=covf, 
                            stderr=STDOUT
                            )

def amptable(ampstats: str):

    odf = pd.DataFrame(columns=[
        'amplicon_number',
        'amplicon_reads',
        'full_length_depth',
        'avg_depth',
        'amplicon_coverage'
        ]
    )
    with open(ampstats, 'r') as sam:
        for line in sam:
            line = line.strip()
            if line.startswith("FREADS"):
                columns = line.split('\t')[2:]
                odf['amplicon_reads'] = columns
                odf['amplicon_number'] = odf.index + 1
            if line.startswith("FVDEPTH"):
                columns = line.split('\t')[2:]
                odf['full_length_depth'] = columns
            if line.startswith("FDEPTH"):
                columns = line.split('\t')[2:]
                odf['avg_depth'] = columns
            if line.startswith("FPCOV-1"):
                columns = line.split('\t')[2:]
                odf['amplicon_coverage'] = columns

    return odf

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

