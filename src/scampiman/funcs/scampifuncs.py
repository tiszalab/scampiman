import logging
import re, shutil
import pysam
import sys, os, time
import mappy as mp
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
    elif job == "bam":
        td = timedelta(seconds=int(time_taken))
        done = elapsed_process == total_process
        mark = " ✔" if done else ""
        sys.stdout.write(f"\r🦐 {percent:3d}% | Processing BAMs{mark} | {td}  ")
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
_sq_head = None
_cpus = 0

def mappy_al_ref(ref: str, tech: str, sq_head: dict, cpus: int):
    """Set up reference and BAM header for the alignment.
    
    Args:
        ref (str): Reference genome file path.
        tech (str): Sequencing technology (ont or illumina).
        sq_head (dict): BAM header dict.
    """
    global _aligner, _sq_head, _ali_count
    if tech == "ont":
        targ = "lr:hq"
    elif tech == "illumina":
        targ = "sr"

    _aligner = mp.Aligner(ref,
        preset=targ)
    if not _aligner:
        raise RuntimeError(f"mappy failed to load/build index from: {ref_path}")
    _sq_head = sq_head
    _cpus = cpus


# Chunking helper
def _chunked(iterable, size: int, outd: str, sample:str) -> Generator[list, None, None]:
    """ Chunking helper 
    Args:
        iterable: Iterable to chunk.
        size (int): Size of each chunk.

    Yields:
        list: Chunk of the iterable.
    """

    chunk: list = []
    chunk_num = 1

    for item in iterable:
        chunk.append(item)
        if len(chunk) == size:
            yield {"chunk": chunk, "chunk_num": chunk_num, "outd": outd, "sample": sample}
            chunk_num += 1
            chunk = []
    if chunk:
        yield {"chunk": chunk, "chunk_num": chunk_num, "outd": outd, "sample": sample}



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
    reads_list = []
    fasta = pysam.FastaFile(ref)
    sq_head = pysam.AlignmentHeader(
    ).from_references(
        reference_names=fasta.references, 
        reference_lengths=fasta.lengths
        ).to_dict()
    fasta.close()
    sq_head['PG'] = []
    if rfmt == "bam":
        for file in file_list:
            bam = pysam.AlignmentFile(file, 'rb', check_sq=False)
            for i in bam:
                reads_list.append((i.query_name, i.query_sequence, i.query_qualities))
                read_count += 1
            if len(sq_head) == 2: # the only thing in the freshly made header is the SQ tag and an empty PG tag
                sq_head = bam.header.to_dict() | sq_head # combine the headers
            else: # if the header has already been combined
                sq_head['RG'] += bam.header.to_dict()['RG'] # add in the RG tag just incase smaples were from different runs
            bam.close()
        # get rid of duplicate RG tags from looping
        g = list({frozenset(d.items()) for d in sq_head['RG']}) # convert the list into a set to remove duplicates
        g = [dict(f) for f in g] # Convert back to dictionaries
        sq_head['RG'] = g # set RG header tag
    elif rfmt == "fastq":
        if file2: # paired end
            for fastq1, fastq2 in zip(file_list, file2): # zip file1 and file2 to be processed in parallel
                with pysam.FastxFile(fastq1) as fq1, pysam.FastxFile(fastq2) as fq2: # open fastq files
                    for read1, read2 in zip(fq1, fq2): # zip read1 and read2 to be processed in parallel
                        read_count += 2
                        reads_list.append((read1.name, read1.sequence, read1.quality, read2.name, read2.sequence, read2.quality))
        else: # single end
            for file in file_list:
                for i in pysam.FastxFile(file): 
                    read_count += 1
                    reads_list.append((i.name, i.sequence, i.quality)) # empty list for tags to keep the single end processing consistant
    # add in the PG tag (gives details about mappy and which read files were used for alignment)
    sq_head['PG'].append({'ID': 'aligner', 'PN': 'mappy', 'VN': mp.__version__, 'DS': 'minimap2 alignment', 'CL': ' '.join(file_list)})
    logger.info(f"Total reads to align: {read_count}")
    shrimp_progress(read_count, 1, 0, "processing") # reporting the read number

    sq_head_obj = pysam.AlignmentHeader.from_dict(sq_head) # make the header object
    return sq_head_obj, read_count, reads_list




# formatting the hits 
def _fmt_se_hits(q: dict) -> Tuple[List[Dict[str, Any]], bool]:
    """Pre-compute CIGAR strings, flags, and tags for single-end hits.
    Pure Python / mappy only — no pysam objects — safe to run inside workers.

    Returns:
        (formatted, is_qcfail) where *formatted* is one spec-dict per hit and
        - "formatted" is a list of hit dictionaries containing the fields needed for a pysam formatted alignedsegment.
        - "is_qcfail" is a boolean indicating if the read should be marked to be filtered out
    """
    hits = q['hits']
    seq  = q['query_sequence']
    qual = q['query_qualities']
    sq_head = q['sq_header']
    primary_hit = hits[0]
    is_qcfail = False
    formatted: List[Dict[str, Any]] = []
    primary_rec = None

    for hit in hits:
        spec = {'is_primary': False, 'is_secondary': False, 'is_supplementary': False}
        rec = pysam.AlignedSegment(header=sq_head)
        rec.is_mapped = True
        rec.query_name = q['query_name']
        rec.reference_name = hit.ctg 
        rec.reference_start = hit.r_st 
        rec.mapping_quality = hit.mapq 

        if hit is primary_hit: # only set these tags for the primary hit 
            spec['is_primary'] = True
            if hit.strand == 1:
                rec.is_forward = True
                rec.query_sequence = seq
                rec.query_qualities = qual
                c_front = hit.q_st # which base of the query sequence is the first base of the alignment; how many bases are clipped from the front
                c_end   = len(seq) - hit.q_en # which base of the query sequence is the last base of the alignment; how many bases are clipped from the end
            else: # reverse strand
                rec.is_reverse = True               
                rec.query_sequence  = mp.revcomp(seq) # bam output for reverse strand is reverse complement
                rec.query_qualities = qual[::-1] # the quality string is also reverse
                c_front = len(seq) - hit.q_en # which base of the query sequence (as reverse) is the first base of the alignment
                c_end   = hit.q_st # which base of the query sequence (as reverse) is the last base of the alignment
            if c_front and c_end: # if there are clippings at both ends
                rec.cigarstring = f"{c_front}S{hit.cigar_str}{c_end}S"
            elif c_front: # if there is a clipping at the front but not the back
                rec.cigarstring = f"{c_front}S{hit.cigar_str}"
            elif c_end: # if there is a clipping at the end but not the front
                rec.cigarstring = f"{hit.cigar_str}{c_end}S"
            else: # if there are no clippings
                rec.cigarstring = hit.cigar_str
            primary_rec = rec
        else:
            prime_q   = range(min(primary_hit.q_st, primary_hit.q_en), max(primary_hit.q_st, primary_hit.q_en)) # range of base positions on the read of the primary read alignment
            hit_q     = range(min(hit.q_en, hit.q_st), max(hit.q_en, hit.q_st)) # range of base positions on the read of the hit
            qlen_olap = len(list(set(prime_q) & set(hit_q))) / len(hit_q) # overlap of the hit with the primary read alignment as a fraction of the hit length
            if qlen_olap < 0.5:
                rec.is_supplementary = True
                spec['is_supplementary'] = True
            else:
                rec.is_secondary = True
                spec['is_secondary'] = True

            if hit.strand == 1:
                rec.is_forward = True
                rec.query_sequence  = seq[hit.q_st:hit.q_en] # report only the aligned sequence bases
                rec.query_qualities = qual[hit.q_st:hit.q_en] # report only the aligned quality positions
                rec.cigarstring = f"{hit.q_st}H{hit.cigar_str}{len(seq)-hit.q_en}H" # building cigar string for secondary alignment; H = hard clipping
            else:
                rec.is_forward = False
                rec.query_sequence  = mp.revcomp(seq[hit.q_st:hit.q_en]) # report only the aligned sequence bases as reverse complement
                rec.query_qualities = qual[hit.q_st:hit.q_en][::-1] # report only the aligned quality positions as reverse
                rec.cigarstring = f"{len(seq)-hit.q_en}H{hit.cigar_str}{hit.q_st}H" # building cigar string for secondary alignment; H = hard clipping

            ref_olap = primary_rec.get_overlap(hit.r_st, hit.r_en) # number of overlapping reference bases between the primary hit and the hit
            hit_refRange = hit.r_en - hit.r_st # reference length of the hit
            rlen_olap_ratio = ref_olap/hit_refRange # ratio of overlapping reference bases of the primary and the hit to the length of the hit; how much of the hit is covered by the primary hit
            if rlen_olap_ratio < 0.5: ## if the overlap is less than 50% of the hit, it is a secondary alignment 
                is_qcfail = True
            elif rec.is_supplementary and primary_hit.strand == hit.strand:  # if the hit is supplementary and on the same strand as the primary hit
                is_qcfail = True

       ## creating tags
        ### de tag: Gap-compressed per-base sequence divergence = matching 
        de_denom   = sum(l if op == 0 else 1 for l, op in hit.cigar) # loop through the array of CIGAR [length, operation]; Add up the number of matched bases and any other operation (insertion, deletion, etc.)
        de         = 1 - (hit.mlen / de_denom)
        rounded_de = round(de, str(round(de, 6)).split('.')[1].count('0') + 6)
        
        tp          = 'P' if hit.is_primary else 'S' # tp tag: primary (P) or [supplementary or secondary] (S)
        strand_char = '+' if hit.strand == 1 else '-' # strand character for SA tag: forward (+) or reverse (-)

        spec['tags'] = [
            ("NM", hit.NM),
            ("nn", hit.blen - hit.mlen - hit.NM, "i"),
            ("tp", tp, "A"),
            ("de", rounded_de, "f"),
            ("MD", hit.MD),
            ("cs", hit.cs),
        ]
        spec['rec'] = rec
        spec['sa_tag_str'] = f"{hit.ctg},{hit.r_st},{strand_char},{rec.cigarstring},{hit.mapq},{hit.NM}"
        formatted.append(spec)

    return formatted, is_qcfail


def _fmt_pe_hits(q: dict) -> Tuple[List[Dict[str, Any]], bool]:
    """Pre-compute CIGAR strings, flags, tags, and mate info for paired-end hits.
    Pure Python / mappy only — no pysam objects — safe to run inside workers.

    Returns:
        (formatted, is_qcfail) where *formatted* is one spec-dict per hit
        (plus an extra unmapped-mate spec if one read has no hits) and
        *is_qcfail* indicates the pair should be sent to the fail BAM.
    """


    sq_head = q['sq_header']    
    hits      = q['hits']
    is_qcfail = False
    formatted:  List[Dict[str, Any]] = []
    hit_r1:     List[Dict[str, Any]] = []
    hit_r1_idx: List[int]            = []
    hit_r2:     List[Dict[str, Any]] = []
    hit_r2_idx: List[int]            = []

    for i, hit in enumerate(hits): # the i is the index (enumerate) of the hit in the hits list
        rec = pysam.AlignedSegment(header=sq_head)
        rec.is_mapped = True
        rec.reference_name = hit.ctg 
        rec.reference_start = hit.r_st 
        rec.mapping_quality = hit.mapq
        rec.is_paired = True

        spec: Dict[str, Any] = {
#            "reference_name":       hit['ctg'],
 #           "reference_start":      hit['r_st'],
  #          "mapping_quality":      hit['mapq'],
   #         "is_paired":            True,
            "is_primary":           False,
            "is_supplementary":     False,
            "is_secondary":         False,
            "is_read1":             False,
            "is_read2":             False,
#            "next_reference_name":  None,
#            "next_reference_start": 0,
#            "template_length":      0,
        }

        if hit.read_num == 1:
            spec['is_read1'] = True
            rec.is_read1 = True
            rec.query_name = f"{q['query_name']}_1"
            q_seq = q['query_sequence']
            q_qual = q['query_qualities']
            hit_r1.append(hit)
            hit_r1_idx.append(i) # adds the index of the hit to the hit_r1_idx list
        else:
            spec['is_read2'] = True
            rec.is_read2 = True
            rec.query_name = f"{q['query_name2']}_2"
            q_seq  = q['query_sequence2']
            q_qual = q['query_qualities2']
            hit_r2.append(hit)
            hit_r2_idx.append(i) # adds the index of the hit to the hit_r2_idx list


        if not hit.is_primary:
            tp = 'S'  # set tp tag to 'S' if the hit is not primary
            if hit.read_num == 1:
                prime_q = set(range(hit_r1[0].q_st,hit_r1[0].q_en)) # range of base positions on the read1 of the primary read1 alignment hits[0]
            else:
                prime_q = set(range(hit_r2[0].q_st,hit_r2[0].q_en)) # range of base positions on the read2 of the primary read2 alignment hits[0]
            hit_q     = set(range(hit.q_st, hit.q_en)) # range of base positions on the read of the hit
            qlen_olap = len(prime_q & hit_q) / len(hit_q) # ratio of common READ bases between the primary hit and the hit relative to the hit length
            
            if qlen_olap < 0.5: 
                rec.is_supplementary = True
                spec["is_supplementary"] = True
            else:
                rec.is_secondary = True     
                spec["is_secondary"] = True
        else:
            tp = 'P'
            spec["is_primary"] = True

        if hit.strand == 1:
            rec.is_forward = True
            rec.mate_is_reverse = True
            rec.query_sequence = q_seq
            rec.query_qualities = q_qual
            c_front = hit.q_st # start position of the cigar string on the read (in forward)
            c_end   = len(q_seq) - hit.q_en # end position of the cigar string on the read (in forward)
        else:
            rec.is_reverse = True
            rec.mate_is_forward = True
            rec.query_sequence = mp.revcomp(q_seq) # reverse complement of the read sequence
            rec.query_qualities = q_qual[::-1] # reverse of the read quality/qscore string
            c_front = len(q_seq)-hit.q_en # start position of the cigar string on the read (in reverse)
            c_end = hit.q_st # end position of the cigar string on the read (in reverse)

        # writing cigar string out 
        if c_front and c_end:
            rec.cigarstring = f"{c_front}S{hit.cigar_str}{c_end}S"
        elif c_front:
            rec.cigarstring = f"{c_front}S{hit.cigar_str}"
        elif c_end:
            rec.cigarstring = f"{hit.cigar_str}{c_end}S"
        else:
            rec.cigarstring = hit.cigar_str

        # de tag: Gap-compressed per-base sequence divergence = matching 
        de_denom   = sum(l if op == 0 else 1 for l, op in hit.cigar) # loop through the array of CIGAR [length, operation]; Add up the number of matched bases and any other operation (insertion, deletion, etc.)
        de         = 1 - (hit.mlen / de_denom) 
        rounded_de = round(de, 4)

        rec.tags = [
            ("NM", hit.NM),
            ("nn", hit.blen - hit.mlen - hit.NM, "i"),
            ("tp", tp, "A"),
            ("de", rounded_de, "f"),
            ("MD", hit.MD),
            ("cs", hit.cs),
        ]
        spec['rec'] = rec
        formatted.append(spec)

    # second pass: set mate info, TLEN
    for ((h1, idx1), (h2, idx2)) in zip(
        zip(hit_r1, hit_r1_idx), zip(hit_r2, hit_r2_idx)
    ):
        s1 = formatted[idx1]['rec'] #set the read1 hit
        s2 = formatted[idx2]['rec'] #set the read2 hit

        s1.next_reference_name  = s2.reference_name #set the next reference name for read1
        s1.next_reference_start = s2.reference_start #set the next reference start for read1
        s1.is_proper_pair = True #set the proper pair flag for read1

        s2.next_reference_name  = s1.reference_name #set the next reference name for read2
        s2.next_reference_start = s1.reference_start #set the next reference start for read2
        s2.is_proper_pair = True #set the proper pair flag for read2

        if s1.is_reverse:  # read1 is reverse
            s1.template_length = s1.next_reference_start - s1.reference_start - h1.mlen - h1.NM # read2 start - read1 start - read1 length (not gaps, ambigious, or mismatches) - read1 NM (mismatches)
            s2.template_length = s1.template_length * -1 # read1 template length * -1
        elif s1.is_forward:  # read1 is forward
            s2.template_length = s2.next_reference_start - s2.reference_start - h2.mlen - h2.NM # read1 start - read2 start - read2 length (not gaps, ambigious, or mismatches) - read2 NM (mismatches)
            s1.template_length = s2.template_length * -1 # read2 template length * -1
    
    ## Start qc fail checks
        if s1.is_forward == s2.is_forward: #if both reads are on the same strand
            is_qcfail = True
        # over lap in reference
        if s1.get_overlap(h2.r_st, h2.r_en) == 0: #if the reference mapping is not overlapping between the two reads
            is_qcfail = True

    if len(hit_r1) != len(hit_r2): #if there are not the same number of hits for both reads
        is_qcfail = True

    # partially unmapped: append a stub spec for the missing read
    if not hit_r1 or not hit_r2:
        is_qcfail = True
        unmapped_read = pysam.AlignedSegment(header=sq_head)
        unmapped_read.is_unmapped = True
        unmapped_read.is_paired = True
        if not hit_r1:
            spec['is_read1'] = True
            unmapped_read.is_read1 = True
            unmapped_read.is_forward = True
            unmapped_read.mate_is_reverse = True
            unmapped_read.query_name = f"{q['query_name']}_R1"
            unmapped_read.query_sequence = q['query_sequence']
            unmapped_read.query_qualities = q['query_qualities']
        else:
            spec['is_read2'] = True
            unmapped_read.is_read2 = True
            unmapped_read.is_reverse = True
            unmapped_read.mate_is_forward = True
            unmapped_read.query_name = f"{q['query_name']}_R2"
            unmapped_read.query_sequence = q['query_sequence2']
            unmapped_read.query_qualities = q['query_qualities2']
        spec['rec'] = rec
        formatted.append(spec)

    return formatted, is_qcfail


# Worker functions (run inside the process pool)
def align_se_chunk(chunk_dict: Dict) -> Dict[str, Any]:
    """Align a list of single-end reads, write records to per-chunk temp BAM files.
    Args:
        chunk_dict (Dict): Dict with keys: chunk, chunk_num, sq_head, outd, sample.
    Returns:
        Dict[str, Any]: Dict with outf, foutf (temp BAM paths) and chunk_stats.
    """
    global _aligner, _sq_head
    chunk    = chunk_dict["chunk"]
    outf     = os.path.join(chunk_dict["outd"], f'{chunk_dict["sample"]}.temp{chunk_dict["chunk_num"]}.bam')
    foutf    = os.path.join(chunk_dict["outd"], f'{chunk_dict["sample"]}.temp{chunk_dict["chunk_num"]}.fail.bam')
    obam     = pysam.AlignmentFile(outf,  'wb', header=_sq_head)
    fbam     = pysam.AlignmentFile(foutf, 'wb', header=_sq_head)
    chunk_stats = {'unmapped': 0, 'removed_reads_primary': 0, 'kept_primary': 0, 'kept_secondary': 0, 'kept_supplementary': 0}

    for name, seq, qual in chunk:
        hits = list(_aligner.map(seq, cs=True, MD=True))

        if not hits:
            rec = pysam.AlignedSegment(header=_sq_head)
            rec.query_name      = name
            rec.query_sequence  = seq
            rec.query_qualities = qual
            rec.is_unmapped     = True
            chunk_stats['unmapped'] += 1
            fbam.write(rec)
            continue

        q = {"query_name": name, "query_sequence": seq, "query_qualities": qual, "hits": hits, "sq_header": _sq_head}
        formatted, is_qcfail = _fmt_se_hits(q)

        if is_qcfail:
            chunk_stats['removed_reads_primary'] += 1

        sa_tags = [fmt['sa_tag_str'] for fmt in formatted]
        for hit_index, fmt in enumerate(formatted):
            rec = fmt['rec']
            read_sa_tag = sa_tags.copy()
            read_sa_tag.pop(hit_index)
            rec.tags = fmt['tags'] + [("SA", ';'.join(read_sa_tag))] if read_sa_tag else fmt['tags']

            if is_qcfail:
                fbam.write(rec)
                continue

            if fmt['is_secondary']:
                chunk_stats['kept_secondary'] += 1
            elif fmt['is_supplementary']:
                chunk_stats['kept_supplementary'] += 1
            else:
                chunk_stats['kept_primary'] += 1
            obam.write(rec)

    obam.close()
    fbam.close()
    return {"outf": outf, "foutf": foutf, "chunk_stats": chunk_stats}

def align_pe_chunk(chunk_dict: Dict) -> Dict[str, Any]:
    """Align a list of paired-end reads, write records to per-chunk temp BAM files.
    Args:
        chunk_dict (Dict): Dict with keys: chunk, chunk_num, sq_head, outd, sample.
    Returns:
        Dict[str, Any]: Dict with outf, foutf (temp BAM paths) and chunk_stats.
    """
    global _aligner, _sq_head
    chunk    = chunk_dict["chunk"]
    outf     = os.path.join(chunk_dict["outd"], f'{chunk_dict["sample"]}.temp{chunk_dict["chunk_num"]}.bam')
    foutf    = os.path.join(chunk_dict["outd"], f'{chunk_dict["sample"]}.temp{chunk_dict["chunk_num"]}.fail.bam')
    obam     = pysam.AlignmentFile(outf,  'wb', header=_sq_head)
    fbam     = pysam.AlignmentFile(foutf, 'wb', header=_sq_head)
    chunk_stats = {'unmapped': 0, 'removed_reads_primary': 0, 'kept_primary': 0, 'kept_secondary': 0, 'kept_supplementary': 0, 'read_1': 0, 'read_2': 0}

    for read1_name, read1_sequence, read1_quality, read2_name, read2_sequence, read2_quality in chunk:
        hits = list(_aligner.map(read1_sequence, seq2=read2_sequence, cs=True, MD=True))

        if not hits:
            for rname, rseq, rqual in [(read1_name, read1_sequence, read1_quality), (read2_name, read2_sequence, read2_quality)]:
                rec = pysam.AlignedSegment(header=_sq_head)
                rec.query_name      = rname
                rec.query_sequence  = rseq
                rec.query_qualities = rqual
                rec.is_unmapped     = True
                rec.is_paired       = True
                if rname == read1_name:
                    rec.is_read1 = True
                else:
                    rec.is_read2 = True
                chunk_stats['unmapped'] += 1
                fbam.write(rec)
            continue

        q = {"query_name": read1_name, "query_sequence": read1_sequence, "query_qualities": read1_quality,
             "query_name2": read2_name, "query_sequence2": read2_sequence, "query_qualities2": read2_quality, "hits": hits, "sq_header": _sq_head}
        formatted, is_qcfail = _fmt_pe_hits(q)

        if is_qcfail:
            chunk_stats['removed_reads_primary'] += 2

        for fmt in formatted:
            rec = fmt['rec']

            if is_qcfail:
                fbam.write(rec)
                continue

            obam.write(rec)

            if fmt['is_secondary']:
                chunk_stats['kept_secondary'] += 1
            elif fmt['is_supplementary']:
                chunk_stats['kept_supplementary'] += 1
            else:
                chunk_stats['kept_primary'] += 1
            if fmt['is_read1']:
                chunk_stats['read_1'] += 1
            elif fmt['is_read2']:
                chunk_stats['read_2'] += 1

    obam.close()
    fbam.close()
    return {"outf": outf, "foutf": foutf, "chunk_stats": chunk_stats}



# Do the mappy alignment
def mappy_al(rcon: str, rfmt: str, cpus: int, tech: str, ref: str, outd: str, sample:str, file1:list, file2:list=None):
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
        flagstats['read_1'] = 0
        flagstats['read_2'] = 0
    sq_head, read_count, read_iter = alignment_preprocessing(ref, rfmt, file1) if rcon == "single-end" else alignment_preprocessing(ref, rfmt, file1, file2) # create header and get totalcount of reads
    outf  = os.path.join(outd, f'{sample}.bam')
    foutf = os.path.join(outd, f'{sample}.failed.bam')

    flagstats['total_reads'] = read_count
    chunk_size = int(read_count/(cpus - 1))
    header_starttime = time.perf_counter() # start timer for progress bar
    
    worker_fn = align_se_chunk if rcon == "single-end" else align_pe_chunk

    # Dispatch chunks to the process pool and write results as they arrive
    with ProcessPoolExecutor(
        max_workers=cpus,
        initializer=mappy_al_ref,
        initargs=(ref, tech, sq_head, cpus),
    ) as executor:
        pool_tally = 0
        shrimp_progress(cpus, pool_tally, 0, "align")

        futures = {
            executor.submit(worker_fn, chunk): idx
            for idx, chunk in enumerate(_chunked(read_iter, chunk_size, outd, sample))
        }


        temp_outfs  = []
        temp_foutfs = []

        for future in as_completed(futures):
            pool_tally += 1
            header_endtime = time.perf_counter()
            time_taken = header_endtime - header_starttime
            shrimp_progress(cpus, pool_tally, time_taken, "align")

            result = future.result()
            temp_outfs.append(result["outf"])
            temp_foutfs.append(result["foutf"])
            for key in result["chunk_stats"]:
                if key in flagstats:
                    flagstats[key] += result["chunk_stats"][key]
    bam_starttime = time.perf_counter() # start timer for progress bar
    shrimp_progress(4, 0, (time.perf_counter() - bam_starttime), "bam")
    pysam.cat("-o", outf,  *temp_outfs)
    shrimp_progress(4, 1, (time.perf_counter() - bam_starttime), "bam")
    pysam.cat("-o", foutf, *temp_foutfs)
    shrimp_progress(4, 2, (time.perf_counter() - bam_starttime), "bam")
    for f in temp_outfs + temp_foutfs:
        os.remove(f)
    shrimp_progress(4, 3, (time.perf_counter() - bam_starttime), "bam")
    pysam.sort("-o", outf.replace('.bam', '.sort.bam'), outf)
    shrimp_progress(4, 4, (time.perf_counter() - bam_starttime), "bam")

    return flagstats

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
    f.close()


def amptable(ampstats: str, sampid: str, outf: str):
    """
    Take in the output of ampliconstats, grab the important parts, save it as a tsv file
    """
    # set up the dictionary
    odf = {'accession': [], 'amplicon_number': [], 'lprimer': [], 'rprimer': []}
    # get amplicon information for empty keys within dictionary
    with open(ampstats, 'r') as lhf:
        for line in lhf:
            line = line.strip()
            if line.startswith("AMPLICON"): # this is set up more like a table rows
                ampinfo = line.split('\t')[1:]
                odf['accession'].append(ampinfo[0])
                odf['amplicon_number'].append(ampinfo[1])
                odf['lprimer'].append(ampinfo[2])
                odf['rprimer'].append(ampinfo[3])

    row_num = len(odf['amplicon_number']) # get the number of amplicons/ number of rows
    # reopen file and pull out the rest of the information (the data is set up more like columns)
    with open(ampstats, 'r') as lhf:
        for line in lhf:
            line = line.strip()

            if line.startswith("FREADS"):
                columns = line.split('\t')[2:]
                if row_num == len(columns):
                    odf['amplicon_reads'] = columns
                else:
                    raise Exception("amplicon_reads data incorrect length")
            if line.startswith("FVDEPTH"):
                columns = line.split('\t')[2:]
                if row_num == len(columns):
                    odf['full_length_depth'] = columns
                else:
                    raise Exception("full_length_depth data incorrect length")
            if line.startswith("FDEPTH"):
                columns = line.split('\t')[2:]
                if row_num == len(columns):
                    odf['avg_depth'] = columns
                else:
                    raise Exception("avg_depth data incorrect length")
            if line.startswith("FPCOV-1"):
                columns = line.split('\t')[2:]
                if row_num == len(columns):
                    odf['amplicon_coverage'] = columns
                else:
                    raise Exception("amplicon_coverage data incorrect length")
    #add in the sample_ID column at the end
    odf['sample_ID'] = [sampid for i in range(row_num)] 
    # make tsv file
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
    return max_length * 2

