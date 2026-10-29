import logging
import re, shutil
import pysam
import sys, os, time
import mappy as mp
from math import ceil
from datetime import timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Dict, Generator, List, Optional, Tuple

# since I'm logging, I need the same logger as in the main script
logger = logging.getLogger("pct_logger")

def shrimp_header(version: str): 
    
    terminal_size = shutil.get_terminal_size()
    console_width = terminal_size.columns

    if console_width >= 80:

        print("\n🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\n"
        +"🌊🌊🌊🌊🌊🌊🌊🌊🦐🦐🦐🌊🦐🦐🌊🦐🦐🌊🦐🦐🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\n"
        +"🌊🌊🌊🌊🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\n"
        +"🌊🌊🦐🦐🌊🦐🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🦐🌊🦐🦐🌊🦐🦐🦐🌊🌊🌊🌊🌊🌊\n"
        +"🌊🦐🌊🌊🌊🦐🌊🌊🦐🦐🦐🌊🦐🌊🌊🌊🦐🌊🦐🦐🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🌊🌊🦐🌊\n"
        +"🌊🌊🦐🌊🌊🦐🌊🌊🦐🌊🦐🌊🦐🌊🌊🌊🦐🌊🦐🌊🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🦐🦐🌊🦐🦐🌊🦐🌊\n"
        +"🌊🌊🌊🦐🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🌊🌊🦐🌊🦐🌊🌊🌊🌊🌊🦐🌊🌊🌊🦐🌊🦐🌊🦐🌊🦐🌊🦐🦐🌊\n"
        +"🌊🦐🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🦐🌊🦐🌊🌊🦐🌊\n"
        +"🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\n"
        +"🌊🌊🌊🌊🌊🌊🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🦐🦐🌊🌊🦐🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊" 
        + f"v{version}" + "🌊\n")

    else:
        print("🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\n"
        +"🌊🌊 scampiman - " + f"v{version}" + " 🌊🌊\n"
        +"🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊🌊\n")

def shrimp_progress(total_process: int, elapsed_process: int, time_taken: float, job: str) -> None:
    """Display a lightweight progress indicator with a shrimp emoji.
    
    Args:
        total_process: Total number of processes to complete.
        elapsed_process: Number of processes completed.
        time_taken: Time taken so far in seconds.
        job: One of "preprocessing", "align", or "amp".
    """
    percent = int(elapsed_process / total_process * 100) if total_process else 0

    if job == "preprocessing":
        print("🌊  Pulling in reads 🌊  Counting reads  🌊")
    elif job == "processing":
        print(f"\n🌊 Processing {total_process} reads with shrimp power! 🌊\n")
    elif job == "align":
        td = timedelta(seconds=int(time_taken))
        done = elapsed_process == total_process
        mark = " ✔" if done else ""
        sys.stdout.write(f"\r🦐 {percent:3d}% | Alignment{mark} | {td}  ")
        if done:
            sys.stdout.write("\n")
        sys.stdout.flush()
    elif job == "amp":
        td = timedelta(seconds=int(time_taken))
        done = elapsed_process == total_process
        mark = " ✔" if done else ""
        sys.stdout.write(f"\r🦐 {percent:3d}% | Amplicon Analysis{mark} | {td}  ")
        if done:
            sys.stdout.write("\n\n\n✨ Shrimp has arrived! Scampiman complete! ✨\n\n")
        sys.stdout.flush() 



def paired_paths(file_list: list):
    """ Return matching R1 and R2 files from a list of read paths

    Args:
        reads: string of space separated paths 
        rcon: should be"paired-end"
    """
    #establish the file patterns
    r1_pattern = re.compile(r'(.+?)(_R1_|_R1\.|.R1\.|_1\.|_1_)(.*)$') # pattern for R1 files
    r2_pattern = re.compile(r'(.+?)(_R2_|_R2\.|.R2\.|_2\.|_2_)(.*)$') # pattern for R2 files
    r1_files = {}
    r2_files = {}
    #sort the files into R1 and R2
    for fpath in file_list:
        fname = os.path.basename(fpath)
        # match the file name to the pattern for the file
        r1_match = r1_pattern.match(fname)
        r2_match = r2_pattern.match(fname)
        # get the file prefix (everything before the R1 or R2) and put into dictionary
        if r1_match:
            sample_key = r1_match.group(1) # file prefix
            r1_files[sample_key] = fpath # put into dictionary
        elif r2_match:
            sample_key = r2_match.group(1) # file prefix
            r2_files[sample_key] = fpath # put into dictionary      
    matched_r1 = []
    matched_r2 = []
    # look for the matching R1 and R2 files by key (file prefix) and add to the matched read1 and read2 list
    e = ''
    for key in r1_files.keys():
        if key in r2_files.keys():
            matched_r1.append(r1_files[key])
            matched_r2.append(r2_files[key])
        else:
            logger.info(f'[ERROR] No R2 pair found for: {r1_files[key]}') # if it doesn't have a matching R2 file, print a message
            e += f'[ERROR] No R2 pair found for {r1_files[key]}\n'

    # look for R2 files that don't have a matching R1 file
    for key in r2_files.keys():
        if key not in r1_files.keys():
            logger.info(f'[ERROR] No R1 pair found for: {r2_files[key]}') # print message  
            e += f'[ERROR] No R1 pair found for {r2_files[key]}\n'

    if not e:
        # return the matched read1 and read2 lists
        logger.info(f"Paired-end files: {matched_r1, matched_r2}")
        return matched_r1, matched_r2
    else:
        raise Exception(f"Read pairs not found.\n{e}")


# Worker-process state (one aligner + header per process, created by the initializer)
_aligner: Optional[mp.Aligner] = None
_rcon = None
_py_head = None
_temp_dir = None

def mappy_al_ref(ref: str, tech: str, rcon: str, sq_head_dict: dict, temp_dir: str):
    """Set up reference, BAM header, and temp dir for worker processes.
    Args:
        ref (str): Reference genome file path.
        tech (str): Sequencing technology (ont or illumina).
        rcon (str): Read configuration (single-end or paired-end).
        sq_head_dict (dict): BAM header dict.
        temp_dir (str): Directory for per-chunk temp BAM files.
    """
    global _aligner, _rcon, _py_head, _temp_dir
    if tech == "ont":
        targ = "lr:hq"
    elif tech == "illumina":
        targ = "sr"

    _aligner = mp.Aligner(ref,
        preset=targ)
    if not _aligner:
        raise RuntimeError(f"mappy failed to load/build index from: {ref}")
    _rcon = rcon
    _py_head = sq_head_dict
    _temp_dir = temp_dir


def alignment_preprocessing(ref: str, rfmt:str, file_list: list, file2: list = None):
    """Takes in the unaligned read files, parses through them to 
        create header, get total read count, and builds a list of lists of read information 
        to iterate through during alignment executions.
    Args:
        ref (str): Reference genome file path.
        rfmt (str): Read format (bam or fastq).
        file_list (list): List of read files.
        file2 (list, optional): List of paired-end read files. Defaults to None.
    Returns:
        tuple: Tuple containing (read_count, reads_list, sq_head)
    """
    shrimp_progress(1, 1, 0, "preprocessing") # Add "Pulling in reads/Counting reads"
    read_count = 0
    max_read_length = 0
    reads_list = []
    fasta = pysam.FastaFile(ref)
    sq_head = pysam.AlignmentHeader(
    ).from_references(
        reference_names=fasta.references, 
        reference_lengths=fasta.lengths
        ).to_dict()
    fasta.close()
    # add in the PG tag (gives details about mappy and which read files were used for alignment)
    sq_head['PG'] = [{'ID': 'aligner', 'PN': 'mappy', 'VN': mp.__version__, 'DS': 'minimap2 alignment', 'CL': ' '.join(file_list)}]
    read_dict = {'name': '', 'flag': '4', 'ref_name': '*',
            'ref_pos': '0', 'map_quality': '0', 'cigar': '*', 
            'next_ref_name': '*', 'next_ref_pos': '0', 'length': '0', 
            'seq': '', 'qual': '', 'tags': []}
    if rfmt == "bam":
        for file in file_list:
            bam = pysam.AlignmentFile(file, 'rb', check_sq=False)
            for i in bam:
                read_count += 1
                read = read_dict.copy()
                read.update({'name': i.query_name, 'seq': i.query_sequence, 'qual': i.query_qualities_str})
                reads_list.append(read)
                if len(i.query_sequence) > max_read_length:
                        max_read_length = len(i.query_sequence)
            bam.close()
    elif rfmt == "fastq":
        if file2: # paired end
            for fastq1, fastq2 in zip(file_list, file2): # zip file1 and file2 to be processed in parallel
                with pysam.FastxFile(fastq1) as fq1, pysam.FastxFile(fastq2) as fq2: # open fastq files
                    for read1, read2 in zip(fq1, fq2): # zip read1 and read2 to be processed in parallel
                        read_count += 2
                        read1_dict = read_dict.copy()
                        read2_dict = read_dict.copy()
                        read1_dict.update({'name': read1.name, 'seq': read1.sequence, 'qual': read1.quality, 'flag': '69'}) # 69 = read1(64),  is paired(1), is unmapped(4)
                        read2_dict.update({'name': read2.name, 'seq': read2.sequence, 'qual': read2.quality, 'flag': '133'}) # 133 = read2(128),  is paired(1), is unmapped(4)
                        reads_list.append([read1_dict,read2_dict])
                        if len(read1.sequence) + len(read2.sequence) > max_read_length:
                            max_read_length = len(read1.sequence) + len(read2.sequence)
        else: # single end
            for file in file_list:
                for i in pysam.FastxFile(file): 
                    read_count += 1 
                    read1_dict = read_dict.copy()
                    read1_dict.update({'name': i.name, 'seq': i.sequence, 'qual': i.quality})
                    reads_list.append(read1_dict)
                    if len(i.sequence) > max_read_length:
                            max_read_length = len(i.sequence)

    logger.info(f"Total reads to align: {read_count}")
    logger.info(f"Longest read length: {max_read_length}")
   
    shrimp_progress(read_count, 1, 0, "processing")
    return sq_head, read_count, reads_list, max_read_length


# Chunking helper
def _chunked(iterable, size: int) -> Generator[list, None, None]:
    """ Chunking helper 
    Args:
        iterable: Iterable to chunk.
        size (int): Size of each chunk. 

    Yields:
        list: Chunk of the iterable.
    """
    chunk: list = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

def _calc_flag(props: list) -> int:
    """ Calculate the SAM/BAM flag value based on the combination of alignment properties. Inputs created/determined in _fmt_hits.
    Args:  
        props: a list containing the alignment properties of the given record 
    
    Returns:
        int: calculated sum of SAM/BAM flag values. Based on: https://broadinstitute.github.io/picard/explain-flags.html
    """

    ### Dictionary of alignment property (key) and associated flag code value. 
    flags_coding = {'read_paired': 1, 
                'read_mapped_in_proper_pair': 2,
                'read_unmapped': 4,
                'mate_unmapped': 8,
                'read_reverse_strand': 16,
                'mate_reverse_strand': 32,
                'first_in_pair': 64,
                'second_in_pair': 128,
                'not_primary_alignment': 256,
                'qc_failed' : 512, # aka 'read_fails_platform/vendor_quality_checks'
                #'read_is_PCR_or_optical_duplicate': 1024
                'supplementary_alignment': 2048
                }
    
    ### Get unique properties ; just in case of duplicate properities
    uniq_props = list(set(props))

    ### add up flag code values of unique alignment properties
    flag_total = 0 # Start 
    for key in uniq_props:
        flag_total += flags_coding[key]

    ### Duplicate sanity check
    if len(props) != len(uniq_props):
        seen = set()
        duplicates = set()
        for p in props:
            if p in seen:
                duplicates.add(p)
            else:
                seen.add(p)
        print("\n Duplicate alignment properties detected:" + f"{list(duplicates)}")
    
    return flag_total


# formatting the hits 
def _fmt_hits(q: dict) -> Tuple[List[Dict[str, Any]], bool]:
    """Pre-compute CIGAR strings, flags, and tags for single-end hits.
    Pure Python / mappy only — no pysam objects — safe to run inside workers.

    Returns:
        (formatted, is_qcfail) where *formatted* is one spec-dict per hit and
        - "formatted" is a list of hit dictionaries containing the fields needed for a pysam formatted alignedsegment.
        - "is_qcfail" is a boolean indicating if the read should be marked to be filtered out
    """
    global _rcon
    hits = q['hits']

    is_qcfail = False
    fail_tag = "QC:Z:" 

    rec_list = []
    stats = {'is_secondary': 0, 'is_supplementary': 0}
    if _rcon == "paired-end":
        hit_r1 = []
        hit_r2 = []
        hit_r1_idx = []
        hit_r2_idx = []
    else:
        sa_tag_str = []
        primary_hit = hits[0]

    for i, hit in enumerate(hits): # the i is the index (enumerate) of the hit in the hits list:
    ### items to reset with each loop
        is_primary = False
        flag_list = []


    ### is it primary:
        if _rcon == "single-end":
            rec_dict = q['rec']
            rec = rec_dict.copy() # copy the blank record dictionary to avoid modifying it (just in case there are multiple hits)
            #rec['flag'] = 0
            is_primary = primary_hit
        else: #paired-end
            flag_list += ['read_paired']
     ## set read number in pair flag
            if hit.read_num == 1:
                rec_dict = q['rec']
                rec = rec_dict.copy() # copy the blank record dictionary to avoid modifying it (just in case there are multiple hits)
                #rec.update({'flag': 65}) #flag=paired(1), read1(64)
                flag_list += ['first_in_pair']
                #is_primary = hit.is_primary
                is_primary = hit.is_primary 

                hit_r1.append(hit)
                hit_r1_idx.append(i) # adds the index of the hit to the hit_r1_idx list
            else: # read2
                rec_dict = q['rec2']
                rec = rec_dict.copy() # copy the blank record dictionary to avoid modifying it (just in case there are multiple hits)
                #rec.update({'flag': 129}) #flag=paired(1), read2(128)
                flag_list += ['second_in_pair']
                #is_primary = hit.is_primary # because the hard clippings do not work the same way
                is_primary = hit.is_primary
                hit_r2.append(hit)
                hit_r2_idx.append(i) # adds the index of the hit to the hit_r2_idx list

        rec.update({'ref_name': hit.ctg, 'ref_pos': str(hit.r_st + 1), 'map_quality': str(hit.mapq)})

    ### orientation and clipping
        if hit.strand == 1:
            c_front = hit.q_st # which base of the query sequence is the first base of the alignment; how many bases are clipped from the front
            c_end   = len(rec_dict['seq']) - hit.q_en # which base of the query sequence is the last base of the alignment; how many bases are clipped from the end
            #if _rcon == "paired-end":
                #rec['flag'] += 32 #rec.mate_is_reverse = True
            #    flag_dict['mate_reverse_strand'] = True
        else: # reverse strand
            #rec['flag'] += 16 # rec.is_reverse = True; works for both single- and paired-end
            flag_list += ['read_reverse_strand']
            rec['seq']  = mp.revcomp(rec_dict['seq']) # bam output for reverse strand is reverse complement
            rec['qual'] = rec_dict['qual'][::-1] # the quality string is also reverse
            c_front = len(rec_dict['seq']) - hit.q_en # which base of the query sequence (as reverse) is the first base of the alignment
            c_end   = hit.q_st # which base of the query sequence (as reverse) is the last base of the alignment

        if c_front and c_end: # if there are clippings at both ends
            rec['cigar'] = f"{c_front}S{hit.cigar_str}{c_end}S"
        elif c_front: # if there is a clipping at the front but not the back
            rec['cigar'] = f"{c_front}S{hit.cigar_str}"
        elif c_end: # if there is a clipping at the end but not the front
            rec['cigar'] = f"{hit.cigar_str}{c_end}S"
        else: # if there are no clippings
            rec['cigar'] = hit.cigar_str

    ### format non-primary hits
        if not is_primary:
         ## modify clipping for non-primary single-end read alignments
            if _rcon == "single-end":
                rec['cigar'] = rec['cigar'].replace("S", "H") # non-primary single-end reads are hard clipped
                if hit.strand == 1:
                    rec['seq']  = rec_dict['seq'][hit.q_st:hit.q_en] # report only the aligned sequence bases
                    rec['qual'] = rec_dict['qual'][hit.q_st:hit.q_en] # report only the aligned quality positions
                else:
                    rec['seq']  = mp.revcomp(rec_dict['seq'][hit.q_st:hit.q_en]) # report only the aligned sequence bases as reverse complement
                    rec['qual'] = rec_dict['qual'][hit.q_st:hit.q_en][::-1] # report only the aligned quality positions as reverse
         ## identify supplementary vs secondary alignment
          # primary hit start and end positions on query read 
                pq_st, pq_en = min(primary_hit.q_st, primary_hit.q_en), max(primary_hit.q_st, primary_hit.q_en) # still under single-end condition; max/min ignores strand direction
            else: # paired-end; primary is always the first read hit
                if hit.read_num == 1:
                    pq_st, pq_en = min(hit_r1[0].q_st, hit_r1[0].q_en), max(hit_r1[0].q_st, hit_r1[0].q_en) 
                else:
                    pq_st, pq_en = min(hit_r2[0].q_st, hit_r2[0].q_en), max(hit_r2[0].q_st, hit_r2[0].q_en) 
          # non-primary hit start and end positions on query read 
            hq_st, hq_en = min(hit.q_st, hit.q_en), max(hit.q_st, hit.q_en)
          # calculate ratio of overlap of mapped bps on the query read between primary and non-primary hits
            qlen_olap = max(0, min(pq_en, hq_en) - max(pq_st, hq_st)) / (hq_en - hq_st)
            if qlen_olap < 0.5: # if the mapping bps on the read overlap less than 50% of the hit, it is a supplementary alignment
                #rec['flag'] += 2048
                is_supplementary = True
                flag_list += ['supplementary_alignment']
                stats['is_supplementary'] += 1
            else: #if it is more than 50% read overlap, it is secondary
                #rec['flag'] += 256
                flag_list += ['not_primary_alignment']
                stats['is_secondary'] += 1
                is_qcfail = True
                fail_tag = f"{fail_tag}Secondary alignment;"
         ## identify technical artifacts within single-end reads
            if _rcon == "single-end": #QC filter
                pr_st, pr_en = min(primary_hit.r_st, primary_hit.r_en), max(primary_hit.r_st, primary_hit.r_en)
                hr_st, hr_en = min(hit.r_st, hit.r_en), max(hit.r_st, hit.r_en)
                rlen_olap_ratio = max(0, min(pr_en, hr_en) - max(pr_st, hr_st)) / (hr_en - hr_st)
                if rlen_olap_ratio < 0.5: ## if the overlapping is less than 50% of the primary genomic region, it is a chimeric read
                    is_qcfail = True
                    fail_tag = f"{fail_tag}Chimeric read: mapping disjunct;"
                elif is_supplementary and primary_hit.strand == hit.strand:  # if the hit is supplementary and on the same strand as the primary hit
                    is_qcfail = True
                    fail_tag = f"{fail_tag}Chimeric read: strand congruence;"
    ### creating tags
     ## de tag: Gap-compressed per-base sequence divergence = matching
        de_denom   = sum(l if op == 0 else 1 for l, op in hit.cigar) # loop through the array of CIGAR [length, operation]; Add up the number of matched bases and any other operation (insertion, deletion, etc.)
        de         = 1 - (hit.mlen / de_denom)
        rounded_de = round(de, str(round(de, 6)).split('.')[1].count('0') + 6)
        tp          = 'P' if hit.is_primary else 'S' # tp tag: primary (P) or [supplementary or secondary] (S)
     ## additional tags
        #rec['flag'] = str(rec['flag'])
        rec['tags'] = [
            f"NM:i:{hit.NM}",
            f"nn:i:{hit.blen - hit.mlen - hit.NM}",
            f"tp:A:{tp}",
            f"de:f:{rounded_de}"
        ]
     ## save current approximate CIGAR string to build a list for 'SA' tag in single-end read alignment
        if _rcon == "single-end":
            cig_dict = {}
            for l, op in hit.cigar:
                cig_dict[op] = cig_dict.get(op, 0) + l

            i_val, d_val = cig_dict.get(1, 0), cig_dict.get(2, 0)
            diff = i_val - d_val
            indel = f"{diff}I" if diff > 0 else (f"{-diff}D" if diff < 0 else "")
            m_val = cig_dict.get(0) if min(i_val, d_val) == 0 else cig_dict.get(0) + min(i_val, d_val)
            sa_front = f"{c_front}S" if c_front else ""
            sa_end = f"{c_end}S" if c_end else ""
            c_string = sa_front + str(m_val) + "M" + indel + sa_end

            strand_char = '+' if hit.strand == 1 else '-' # strand character for SA tag: forward (+) or reverse (-)
            sa_tag_str.append(f"{hit.ctg},{rec['ref_pos']},{strand_char},{c_string},{hit.mapq},{hit.NM}")
    ### save flag_list in record 'flag' to calculate proper flag outside of loop (to standardize QC fail across all hits of a given read)
        rec['flag'] = flag_list
    ### adding record to rec_list
        rec_list.append(rec)

### second pass outside of loop:
    ### update single-end records
    if _rcon == "single-end":
        for hit_index, segment in enumerate(rec_list):
         ## add 'SA' tag - list of other supplementary alignments (with approximate CIGAR strings)
            read_sa_tag = [s for j, s in enumerate(sa_tag_str) if j != hit_index]
            if read_sa_tag:
                segment['tags']+= [f"SA:Z:{';'.join(read_sa_tag)}"] 
        ## Standardize QC fail result across all hits
            if is_qcfail:
            # add in fail tag
                segment['tags'] = [fail_tag[:-1]] + segment['tags'] # [:-1] removes the last ';'
            # update flag list with QC_fail
                segment['flag'] += ['qc_failed']
        ## replace flag list with proper formatted value
            segment['flag'] = str(_calc_flag(segment['flag']))
    ### update paired-end records            
    else: 
     ## check for mapping issues, if there are not the same number of hits for both reads
        if len(hit_r1) != len(hit_r2): 
            fail_tag = f"{fail_tag}Mapping inconsitency: different number of hits for read1/read2,"
        ## Update the flags and tags of rec_list
            for segment in rec_list:
                segment['flag'] += ['qc_failed']
            ## if the mate is unmapped update flags and tags
                if not hit_r1 or not hit_r2:
                    segment['flag'] += ['mate_unmapped']
                    segment['tags'] = [f"{fail_tag} mate unmapped"] + segment['tags'] 
                else:
                    segment['tags'] = [fail_tag[:-1]] + segment['tags'] # [:-1] removes the last ';'
            ## replace flag list with proper formatted value
                segment['flag'] = str(_calc_flag(segment['flag']))
        ## add unmapped reads to rec_list with new flags and tags
            if not hit_r1:
                rec = q['rec'].copy()
                rec['tags'] = [f"{fail_tag} only mate (read2) mapped"]
                rec_list.append(rec)
            else:
                rec = q['rec2'].copy()
                rec['tags'] = [f"{fail_tag} only mate (read1) mapped"]
                rec_list.append(rec)
     ## format reads mapped in proper pairs
        else: # same number of hits for both reads
            for ((h1, idx1), (h2, idx2)) in zip(zip(hit_r1, hit_r1_idx), zip(hit_r2, hit_r2_idx)): # link pairs
            ## symlink appropriate records in rec_list to modify
                s1 = rec_list[idx1] #set the read1 hit
                s2 = rec_list[idx2] #set the read2 hit
            ## set proper pair flag
                s1['flag'] += ['read_mapped_in_proper_pair']
                s2['flag'] += ['read_mapped_in_proper_pair']
            ### set mate info
                s1['next_ref_name']  = s2['ref_name'] #set the next reference name for read1
                s1['next_ref_pos'] = s2['ref_pos'] #set the next reference start for read1
                s2['next_ref_name']  = s1['ref_name'] #set the next reference name for read2
                s2['next_ref_pos'] = s1['ref_pos'] #set the next reference start for read2
            ### set TLEN and strand information (flag and positions)
                f_start = None
                f_end = None
                r_start = None
                r_end = None
                if h1.strand != 1:  # read1 is reverse
                    s1['length'] = int(s1['next_ref_pos']) - int(s1['ref_pos']) - h1.mlen - h1.NM # read2 start - read1 start - read1 length (not gaps, ambigious, or mismatches) - read1 NM (mismatches)
                    s2['length'] = str(s1['length'] * -1) # read1 template length * -1
                    s1['length'] = str(s1['length'])
                    if h2.strand == 1:
                        s2['flag'] += ['mate_reverse_strand']
                    # for qc check
                    f_start = h2.r_st # for qc_check
                    f_end = h2.r_en # for qc_check
                    r_start = h1.r_st # for qc_check
                    r_end = h1.r_en # for qc_check

                elif h1.strand == 1:  # read1 is forward
                    s2['length']= int(s2['next_ref_pos']) - int(s2['ref_pos']) - h2.mlen - h2.NM # read1 start - read2 start - read2 length (not gaps, ambigious, or mismatches) - read2 NM (mismatches)
                    s1['length'] = str(s2['length'] * -1) # read2 template length * -1
                    s2['length'] = str(s2['length'])
                    if h2.strand != 1:
                        s1['flag'] += ['mate_reverse_strand']
                    # for qc check
                    f_start = h1.r_st
                    f_end = h1.r_en
                    r_start = h2.r_st
                    r_end = h2.r_en

            ### assumming pass QC filter, set flag value
                s1['flag'] = str(_calc_flag(s1['flag']))
                s2['flag'] = str(_calc_flag(s2['flag']))
            ## Start paired qc fail checks and add tags
                paired_ftags = fail_tag
                if h1.strand == h2.strand: #if both reads are on the same strand
                    is_qcfail = True
                    paired_ftags = f"{paired_ftags}Same strand orientation of read1/read2;"

                # over lap in reference
                ref_olap = f_end - r_start
                #if the reference mapping is not overlapping between the two reads
                if ref_olap < 0: 
                    is_qcfail = True
                    paired_ftags = f"{paired_ftags}No mapped overlap between read1/read2;"

                if f_start >= r_start or f_end >= r_end:
                    is_qcfail = True
                    paired_ftags = f"{paired_ftags}Primer dimer;"

            ### Set paired failed tags 
                if is_qcfail:
                    s1['tags'] = [paired_ftags[:-1]] + s1['tags']
                    s2['tags'] = [paired_ftags[:-1]] + s2['tags']
        ## Standardize QC fail result across all hits
            if is_qcfail:
                for segment in rec_list:
                    segment['flag'] = str(int(segment['flag']) + 512) # 512 = qc_failed code


    return {'recs': rec_list, 'rec_stats': stats}, is_qcfail


# Worker functions (run inside the process pool)
def align_chunk(chunk: list, chunk_idx: int) -> Dict[str, Any]:
    """Align a list of reads, write results to per-chunk temp BAM files.
    Args:
        chunk (list): Read dicts (SE) or [read1_dict, read2_dict] pairs (PE).
        chunk_idx (int): Unique index for naming temp BAM files.
    Returns:
        Dict[str, Any]: pass_bam, fail_bam paths and chunk_stats.
    """
    global _aligner, _rcon, _py_head, _temp_dir
    chunk_stats = {'unmapped':0,'removed_reads_primary':0, 'kept_primary':0, 'kept_secondary':0, 'kept_supplementary':0} # create flagstats dictionary
    if _rcon == "paired-end":
        chunk_stats.update({'read1':0, 'read2':0})

    py_head = pysam.AlignmentHeader.from_dict(_py_head)
    bam_path = os.path.join(_temp_dir, f"chunk_{chunk_idx}.bam")
    obam = pysam.AlignmentFile(bam_path, 'wb', header=py_head)

    for read_info in chunk:
        if _rcon == "single-end":
            hits = list(_aligner.map(read_info['seq']))
            q = {'rec': read_info, "hits": hits}
            reads = [read_info]
        else:
            read1_dict = read_info[0]
            read2_dict = read_info[1]
            hits = list(_aligner.map(read1_dict['seq'], seq2=read2_dict['seq']))
            q = {'rec': read1_dict, 'rec2': read2_dict, 'hits': hits}
            reads = [read1_dict, read2_dict]

        if not hits:
            for i in reads:
                chunk_stats['unmapped'] += 1
                obam.write(pysam.AlignedSegment.from_dict(i, header=py_head))
            continue

        formatted, is_qcfail = _fmt_hits(q)

        recs = formatted['recs']
        specs = formatted['rec_stats']

        for rec in recs:
            obam.write(pysam.AlignedSegment.from_dict(rec, header=py_head))

        if is_qcfail:
            chunk_stats['removed_reads_primary'] += 1 if _rcon == "single-end" else 2
        else:
            chunk_stats['kept_primary'] += 1 if _rcon == "single-end" else 2
            chunk_stats['kept_secondary'] += specs['is_secondary']
            chunk_stats['kept_supplementary'] += specs['is_supplementary']
            if _rcon == "paired-end":
                chunk_stats['read1'] += 1
                chunk_stats['read2'] += 1

    obam.close()
    return {'bam_path': bam_path, 'chunk_stats': chunk_stats}


# Do the mappy alignment
def mappy_al(rcon: str, rfmt: str, cpus: int, tech: str, amp_length: int, ref: str, outf: str, file1:list, file2:list=None):
    """Align reads to reference.
        Args:
        rcon: read connection type (single-end/paired-end)
        rfmt: read format (bam/fastq)
        cpus: number of cpus
        tech: sequencing technology (ont/illumina)
        ref: reference file
        outf: output file
        foutf: output file for failed reads
        file1: list of files to align
        Returns: alignment stats as pandas dataframe
    """
    flagstats = {'total_reads':0, 'unmapped':0,'removed_reads_primary':0, 'kept_primary':0, 'kept_secondary':0, 'kept_supplementary':0} # create flagstats dictionary
    if rcon == "paired-end":
        flagstats.update({'read1':0, 'read2':0})

    sq_head, read_count, read_iter, read_length = alignment_preprocessing(ref, rfmt, file1) if rcon == "single-end" else alignment_preprocessing(ref, rfmt, file1, file2) # create header and get totalcount of reads
    flagstats['total_reads'] = read_count

    e = ''
    if read_length < amp_length:
        logger.info(f'[ERROR] Read length insufficent for amplicon size: Max read length {read_length} bp') # print message  
        e += f'[ERROR] Read length insufficent for amplicon size: Max read length {read_length} bp \n'

    if not e:

        temp_dir = os.path.dirname(outf)

        num_items = len(read_iter) # number of items to chunk (SE: reads, PE: read pairs)
        chunk_size = 10000 if max(1, ceil(num_items/10000)) < 150 else ceil(num_items/150)
        num_chunks = max(1, ceil(num_items/chunk_size))

        align_starttime = time.perf_counter() # start timer for progress bar

        bam_fpaths = []

        # Dispatch chunks to the process pool; workers write temp BAMs directly
        with ProcessPoolExecutor(
            max_workers=cpus,
            initializer=mappy_al_ref,
            initargs=(ref, tech, rcon, sq_head, temp_dir),
        ) as executor:
            pool_tally = 0
            shrimp_progress(num_chunks, pool_tally, 0, "align")

            futures = {
                executor.submit(align_chunk, chunk, idx): idx
                for idx, chunk in enumerate(_chunked(read_iter, chunk_size))
            }

            for future in as_completed(futures):
                result = future.result()
                bam_fpaths.append(result['bam_path'])

                for key in result["chunk_stats"]:
                    if key in flagstats:
                        flagstats[key] += result["chunk_stats"][key]

                pool_tally += 1
                align_endtime = time.perf_counter()
                time_taken = align_endtime - align_starttime
                shrimp_progress(num_chunks, pool_tally, time_taken, "align")

        # Concatenate per-chunk BAMs into final output files
        pysam.cat("-o", outf, *bam_fpaths)
        pysam.samtools.sort("-@", str(cpus), "-o", outf.replace('.bam', '.sort.bam'), outf)
        # Clean up temp chunk BAMs
        for f in bam_fpaths:
            os.remove(f)

        return flagstats

    else:
        raise Exception(f"{e}")


def process_ampclip_bams(bam_path: str, failed_path: str, cpus: int):
    """
    Process amplicon clip BAM files to identify reads that failed alignment.
    """
    remove_path = os.path.join(os.path.dirname(failed_path), f"remove_unclipped_reads.txt")
    check = set()
    with open(remove_path, 'w') as out:
        for reads in pysam.samtools.view("-F", "512", failed_path, catch_stdout=True).strip().split('\n'):
            name = reads.strip().split('\t')[0]
            if name in check:
                continue
            else:
                out.write(name + '\n')
                check.add(name)
    
    pysam.samtools.view("-@", str(cpus), "-bS", "-N", f"^{remove_path}", "-o", bam_path.replace('.temp.ampclip.bam', '.temp2.ampclip.bam'), bam_path, catch_stdout=False)
    pysam.samtools.sort("-@", str(cpus), "-o", bam_path.replace('.temp.ampclip.bam', '.ampclip.bam'), bam_path.replace('.temp.ampclip.bam', '.temp2.ampclip.bam'))
    
    os.remove(bam_path.replace('.temp.ampclip.bam', '.temp2.ampclip.bam'))
    os.remove(bam_path)


def stats_tsv(stat_dict: dict, nrows: int, outf: str):
    """
    Write stats in a dictionary into a tsv file
    """
    # set up rows for tsv file
    col_names = []
    rows = [[] for i in range(nrows)] # a list of empty list (1 list per row)

    # got through the dictionary and buid the file rows
    for key, value in stat_dict.items(): # makes sure they are in the correct order
        col_names.append(key) # get the column names
        # go through the values of value items and append it to the proper row list
        ## this is basically transposing the lists (if there is one)
        if type(value) is list:
            for r in range(len(value)):
                rows[r].append(value[r]) 
        elif type(value) is int: # if it is just a single value, there is only 1 row. 
            rows[0].append(str(value)) # append it to the first row

    with open(outf, 'w') as f:
        f.write('\t'.join(col_names) + '\n')
        for wr in rows: # write the rows
            f.write('\t'.join(wr) + '\n')


def amptable(ampstats: str, sampid: str, outf: str):
    """
    Take in the output of ampliconstats, grab the important parts, save it as a tsv file
    """
    odf = {'accession': [], 'amplicon_number': [], 'lprimer': [], 'rprimer': []}
    ampstats_cols = {
        'FREADS': 'amplicon_reads',
        'FVDEPTH': 'full_length_depth',
        'FDEPTH': 'avg_depth',
        'FPCOV-1': 'amplicon_coverage',
    }
    with open(ampstats, 'r') as lhf:
        for line in lhf:
            fields = line.rstrip('\n').split('\t')
            if fields[0] == 'AMPLICON':
                odf['accession'].append(fields[1])
                odf['amplicon_number'].append(fields[2])
                odf['lprimer'].append(fields[3])
                odf['rprimer'].append(fields[4])
            elif fields[0] in ampstats_cols:
                odf[ampstats_cols[fields[0]]] = fields[2:]
 
    row_num = len(odf['amplicon_number'])
    for col, name in ampstats_cols.items():
        if name not in odf or len(odf[name]) != row_num:
            raise Exception(f"{name} data incorrect length")
 
    odf['sample_ID'] = [sampid] * row_num
    stats_tsv(odf, row_num, outf)



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


def amplicon_length(bed: str) -> int:
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
    max_length = -1
    for coords in amplicons.values():
        if 'left_start' not in coords or 'right_end' not in coords:
            continue
        length = coords['right_end'] - coords['left_start']
        if length > max_length:
            max_length = length
    if max_length == -1:
        raise ValueError('No complete amplicon pairs found in BED file')
    return max_length