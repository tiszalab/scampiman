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


# Worker-process state (one aligner per process, created by the initializer)
_aligner: Optional[mp.Aligner] = None

def mappy_al_ref(ref: str, tech: str):
    """Set up reference for the alignment."""
    global _aligner
    if tech == "ont":
        targ = "lr:hq"
    elif tech == "illumina":
        targ = "sr"

    _aligner = mp.Aligner(ref,
        preset=targ)
    if not _aligner:
        raise RuntimeError(f"mappy failed to load/build index from: {ref_path}")


# get header, read count and reads iteratored prior to alignment
def alignment_preprocessing(ref: str, rfmt:str, file_list: list, file2: list = None):
    """Make a header for the alignment file."""
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
                reads_list.append((i.query_name, i.query_sequence, i.query_qualities, i.get_tags()))
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
        if file2:
            for fastq1, fastq2 in zip(file_list, file2): # zip file1 and file2 to be processed in parallel
                with pysam.FastxFile(fastq1) as fq1, pysam.FastxFile(fastq2) as fq2: # open fastq files
                    for read1, read2 in zip(fq1, fq2): # zip read1 and read2 to be processed in parallel
                        read_count += 2
                        reads_list.append((read1.name, read1.sequence, read1.quality, read2.name, read2.sequence, read2.quality))
                        #reads_list.append((read2.name, read2.sequence, read2.quality, []))
        else:
            for file in file_list:
                for i in pysam.FastxFile(file):
                    read_count += 1
                    reads_list.append((i.name, i.sequence, i.quality, []))
    # add in the PG tag (gives details about mappy and which read files were used for alignment)
    sq_head['PG'].append({'ID': 'aligner', 'PN': 'mappy', 'VN': mp.__version__, 'DS': 'minimap2 alignment', 'CL': ' '.join(file_list)})
    logger.info(f"Total reads to align: {read_count}")
    shrimp_progress(read_count, 1, 0, "processing") # reporting the read number

    sq_head_obj = pysam.AlignmentHeader.from_dict(sq_head)
    return sq_head_obj, read_count, reads_list


# Helpers: convert C-level mappy.Alignment to a serialisable dict
def hit_to_dict(hit: mp.Alignment) -> Dict[str, Any]:
    hit_dict = {
        "ctg":          hit.ctg,
        "ctg_len":      hit.ctg_len,
        "r_st":         hit.r_st,
        "r_en":         hit.r_en,
        "q_st":         hit.q_st,
        "q_en":         hit.q_en,
        "strand":       hit.strand,
        "mapq":         hit.mapq,
        "blen":         hit.blen,
        "mlen":         hit.mlen,
        "NM":           hit.NM,
        "is_primary":   hit.is_primary,
        "cigar_str":    hit.cigar_str,
        "cigar":        hit.cigar,
        "cs":           hit.cs,
        "MD":           hit.MD,
        "trans_strand": getattr(hit, "trans_strand", None),
    }
    if hit.read_num:
        hit_dict["read_num"] = hit.read_num
        #hit_dict["is.paired"] = hit.proper_pair
    return hit_dict


# Chunking helper
def _chunked(iterable, size: int) -> Generator[list, None, None]:
    chunk: list = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

# Worker functions (run inside the process pool)
def align_se_chunk(chunk: List) -> List[Dict[str, Any]]:
    """Align a list of single-end reads and return serialisable results."""
    global _aligner
    out: List[Dict[str, Any]] = []
    for name, seq, qual, original_tags in chunk:
        hits = [hit_to_dict(h) for h in _aligner.map(seq, cs=True, MD=True)]
        q = {"query_name": name, "query_sequence": seq, "query_qualities": qual, "hits": hits}
        if hits:
            formatted, is_qcfail = _fmt_se_hits(q)
        else:
            formatted, is_qcfail = [], False
        out.append({
            "query_name":      name,
            "query_sequence":  seq,
            "query_qualities": qual,
            "original_tags":   original_tags,
            "hits":            hits,
            "formatted":       formatted,
            "is_qcfail":       is_qcfail,
        })
    return out

def align_pe_chunk(chunk: List) -> List[Dict[str, Any]]:
    """Align a list of single-end reads and return serialisable results."""
    global _aligner
    out: List[Dict[str, Any]] = []
    for read1_name, read1_sequence, read1_quality, read2_name, read2_sequence, read2_quality in chunk:
        hits = [hit_to_dict(h) for h in _aligner.map(read1_sequence, seq2=read2_sequence, cs=True, MD=True)]
        q = {"query_name": read1_name, "query_sequence": read1_sequence, "query_qualities": read1_quality, "query_name2": read2_name, "query_sequence2": read2_sequence, "query_qualities2": read2_quality, "hits": hits}
        if hits:
            formatted, is_qcfail = _fmt_pe_hits(q)
        else:
            formatted, is_qcfail = [], False
        out.append({
            "query_name":      read1_name,
            "query_sequence":  read1_sequence,
            "query_qualities": read1_quality,
            "query_name2":      read2_name,
            "query_sequence2":  read2_sequence,
            "query_qualities2": read2_quality,
            "original_tags":   [],
            "hits":            hits,
            "formatted":       formatted,
            "is_qcfail":       is_qcfail,
        })
    return out


# formatting the hits
def _fmt_se_hits(q: dict) -> Tuple[List[Dict[str, Any]], bool]:
    """Pre-compute CIGAR strings, flags, and tags for single-end hits.
    Pure Python / mappy only — no pysam objects — safe to run inside workers.

    Returns:
        (formatted, is_qcfail) where *formatted* is one spec-dict per hit and
        *is_qcfail* indicates the primary alignment should be sent to the fail BAM.
    """
    hits = q['hits']
    seq  = q['query_sequence']
    qual = q['query_qualities']
    primary_hit = hits[0]
    is_qcfail = False
    formatted: List[Dict[str, Any]] = []

    for hit in hits:
        spec: Dict[str, Any] = {
            "query_name": q['query_name'],
            "reference_name":  hit['ctg'],
            "reference_start": hit['r_st'],
            "mapping_quality": hit['mapq'],
            "is_supplementary": False,
            "is_secondary":     False,
        }

        if hit is primary_hit: # only set these tags for the primary hit 
            if hit['strand'] == 1:
                spec['is_forward'] = True
                spec['is_reverse'] = False
                spec['seq']  = seq
                spec['qual'] = qual
                c_front = hit['q_st'] # which base of the query sequence is the first base of the alignment; how many bases are clipped from the front
                c_end   = len(seq) - hit['q_en'] # which base of the query sequence is the last base of the alignment; how many bases are clipped from the end
            else: # reverse strand
                spec['is_forward'] = False
                spec['is_reverse'] = True               
                spec['seq']  = mp.revcomp(seq) # bam output for reverse strand is reverse complement
                spec['qual'] = qual[::-1] # the quality string is also reverse
                c_front = len(seq) - hit['q_en'] # which base of the query sequence (as reverse) is the first base of the alignment
                c_end   = hit['q_st'] # which base of the query sequence (as reverse) is the last base of the alignment
            if c_front and c_end: # if there are clippings at both ends
                spec['cigarstring'] = f"{c_front}S{hit['cigar_str']}{c_end}S"
            elif c_front:
                spec['cigarstring'] = f"{c_front}S{hit['cigar_str']}"
            elif c_end:
                spec['cigarstring'] = f"{hit['cigar_str']}{c_end}S"
            else:
                spec['cigarstring'] = hit['cigar_str']

        else:
            prime_q   = range(min(primary_hit['q_st'], primary_hit['q_en']), max(primary_hit['q_st'], primary_hit['q_en']))
            hit_q     = range(min(hit['q_en'], hit['q_st']), max(hit['q_en'], hit['q_st']))
            qlen_olap = len(list(set(prime_q) & set(hit_q))) / len(hit_q) #if hit_q else 0
            if qlen_olap < 0.5:
                spec['is_supplementary'] = True
            else:
                spec['is_secondary'] = True

            if hit['strand'] == 1:
                spec['is_forward'] = True
                spec['is_reverse'] = False
                spec['seq']  = seq[hit['q_st']:hit['q_en']]
                spec['qual'] = qual[hit['q_st']:hit['q_en']]
                spec['cigarstring'] = f"{hit['q_st']}H{hit['cigar_str']}{len(seq)-hit['q_en']}H"
            else:
                spec['is_reverse'] = True
                spec['is_forward'] = False
                spec['seq']  = mp.revcomp(seq[hit['q_st']:hit['q_en']])
                spec['qual'] = qual[hit['q_st']:hit['q_en']][::-1]
                spec['cigarstring'] = f"{len(seq)-hit['q_en']}H{hit['cigar_str']}{hit['q_st']}H"

            prim_refRange      = range(min(primary_hit['r_st'], primary_hit['r_en']), max(primary_hit['r_st'], primary_hit['r_en']))
            hit_refRange  = range(min(hit['r_en'], hit['r_st']), max(hit['r_en'], hit['r_st']))
            rlen_olap_ratio = len(list(set(prim_refRange) & set(hit_refRange))) / len(hit_refRange) #if hit_refRange else 0
            if rlen_olap_ratio < 0.5:
                is_qcfail = True
            elif spec['is_supplementary'] and primary_hit['strand'] == hit['strand']:
                is_qcfail = True

        de_denom   = sum(l if op == 0 else 1 for l, op in hit['cigar'])
        de         = 1 - (hit['mlen'] / de_denom)
        rounded_de = round(de, str(round(de, 6)).split('.')[1].count('0') + 6)
        tp          = 'P' if hit['is_primary'] else 'S'
        strand_char = '+' if hit['strand'] == 1 else '-'

        spec['tags'] = [
            ("NM", hit['NM']),
            ("nn", hit['blen'] - hit['mlen'] - hit['NM'], "i"),
            ("tp", tp, "A"),
            ("de", rounded_de, "f"),
            ("MD", hit['MD']),
            ("cs", hit['cs']),
        ]
        spec['sa_tag_str'] = f"{hit['ctg']},{hit['r_st']},{strand_char},{spec['cigarstring']},{hit['mapq']},{hit['NM']}"
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
    hits      = q['hits']
    is_qcfail = False
    formatted:  List[Dict[str, Any]] = []
    hit_r1:     List[Dict[str, Any]] = []
    hit_r1_idx: List[int]            = []
    hit_r2:     List[Dict[str, Any]] = []
    hit_r2_idx: List[int]            = []

    for i, hit in enumerate(hits): # the i is the index (enumerate) of the hit in the hits list
        spec: Dict[str, Any] = {
            "reference_name":       hit['ctg'],
            "reference_start":      hit['r_st'],
            "mapping_quality":      hit['mapq'],
            "is_paired":            True,
            "is_supplementary":     False,
            "is_secondary":         False,
            "is_proper_pair":       False,
            "is_unmapped":          False,
            "next_reference_name":  None,
            "next_reference_start": 0,
            "template_length":      0,
        }

        if hit['read_num'] == 1:
            spec['is_read1']   = True
            spec['is_read2']   = False
            spec['query_name'] = f"{q['query_name']}_1"
            q_seq  = q['query_sequence']
            q_qual = q['query_qualities']
            hit_r1.append(hit)
            hit_r1_idx.append(i) # adds the index of the hit to the hit_r1_idx list
        else:
            spec['is_read1']   = False
            spec['is_read2']   = True
            spec['query_name'] = f"{q['query_name2']}_2"
            q_seq  = q['query_sequence2']
            q_qual = q['query_qualities2']
            hit_r2.append(hit)
            hit_r2_idx.append(i) # adds the index of the hit to the hit_r2_idx list

        if not hit['is_primary']:
            tp = 'S'
            if hit['read_num'] == 1:
                prime_q = set(range(hit_r1[0]['q_st'], hit_r1[0]['q_en']))
            else: #hit['read_num'] == 2
                prime_q = set(range(hit_r2[0]['q_st'], hit_r2[0]['q_en']))
            hit_q     = set(range(hit['q_st'], hit['q_en']))
            qlen_olap = len(prime_q & hit_q) / len(hit_q)
            if qlen_olap < 0.5:
                spec['is_supplementary'] = True
            else:
                spec['is_secondary'] = True
        else:
            tp = 'P'

        if hit['strand'] == 1:
            spec['is_forward']      = True
            spec['is_reverse']      = False
            spec['mate_is_reverse'] = True
            spec['mate_is_forward'] = False
            spec['seq']  = q_seq
            spec['qual'] = q_qual
            c_front = hit['q_st']
            c_end   = len(q_seq) - hit['q_en']
        else:
            spec['is_forward']      = False
            spec['is_reverse']      = True
            spec['mate_is_reverse'] = False
            spec['mate_is_forward'] = True
            spec['seq']  = mp.revcomp(q_seq)
            spec['qual'] = q_qual[::-1]
            c_front = len(q_seq) - hit['q_en']
            c_end   = hit['q_st']

        if c_front and c_end:
            spec['cigarstring'] = f"{c_front}S{hit['cigar_str']}{c_end}S"
        elif c_front:
            spec['cigarstring'] = f"{c_front}S{hit['cigar_str']}"
        elif c_end:
            spec['cigarstring'] = f"{hit['cigar_str']}{c_end}S"
        else:
            spec['cigarstring'] = hit['cigar_str']

        de_denom   = sum(l if op == 0 else 1 for l, op in hit['cigar'])
        de         = 1 - (hit['mlen'] / de_denom) #if de_denom else 0
        rounded_de = round(de, 4)
        spec['tags'] = [
            ("NM", hit['NM']),
            ("nn", hit['blen'] - hit['mlen'] - hit['NM'], "i"),
            ("tp", tp, "A"),
            ("de", rounded_de, "f"),
            ("MD", hit['MD']),
            ("cs", hit['cs']),
        ]
        formatted.append(spec)

    # second pass: set mate info, TLEN
    for ((h1, idx1), (h2, idx2)) in zip(
        zip(hit_r1, hit_r1_idx), zip(hit_r2, hit_r2_idx)
    ):
        s1 = formatted[idx1]
        s2 = formatted[idx2]

        s1['next_reference_name']  = s2['reference_name']
        s1['next_reference_start'] = s2['reference_start']
        s1['is_proper_pair'] = True

        s2['next_reference_name']  = s1['reference_name']
        s2['next_reference_start'] = s1['reference_start']
        s2['is_proper_pair'] = True

        if s1['is_reverse']:  # read1 is reverse
            s1['template_length'] = s1['next_reference_start'] - s1['reference_start'] - h1['mlen'] - h1['NM']
            s2['template_length'] = s1['template_length'] * -1
        elif s1['is_forward']:  # read1 is forward
            s2['template_length'] = s2['next_reference_start'] - s2['reference_start'] - h2['mlen'] - h2['NM']
            s1['template_length'] = s2['template_length'] * -1

        if s1['is_forward'] == s2['is_forward']:  # same strand
            is_qcfail = True

        s1_refRange = range(min(h1['r_st'], h1['r_en']), max(h1['r_st'], h1['r_en']))
        s2_refRange = range(min(h2['r_en'], h2['r_st']), max(h2['r_en'], h2['r_st']))
        ref_olap = len(list(set(s1_refRange) & set(s2_refRange))) / len(s2_refRange)
        if ref_olap == 0:  # no reference overlap
            is_qcfail = True

    if len(hit_r1) != len(hit_r2): #if there are not the same number of hits for both reads
        is_qcfail = True

    # partially unmapped: append a stub spec for the missing read
    if not hit_r1 or not hit_r2:
        is_qcfail = True
        if not hit_r1:
            unmapped_spec = {
                "query_name": f"{q['query_name']}_R1",
                "seq":        q['query_sequence'],
                "qual":       q['query_qualities'],
                "is_read1":   True,
                "is_read2":   False,
                "is_forward":           True,
                "is_reverse":           False,
                "mate_is_forward":      False,
                "mate_is_reverse":      True,
            }
        else:
            unmapped_spec = {
                "query_name": f"{q['query_name2']}_R2",
                "seq":        q['query_sequence2'],
                "qual":       q['query_qualities2'],
                "is_read1":   False,
                "is_read2":   True,
                "is_forward":           False,
                "is_reverse":           True,
                "mate_is_forward":      True,
                "mate_is_reverse":      False,
            }
        unmapped_spec.update({
            "reference_name":       None,
            "reference_start":      0,
            "mapping_quality":      0,
            "is_paired":            True,
            "is_supplementary":     False,
            "is_secondary":         False,
            "is_proper_pair":       False,
            "is_unmapped":          True,
            "next_reference_name":  None,
            "next_reference_start": 0,
            "template_length":      0,
            "cigarstring":          None,
            "tags":                 [],
        })
        formatted.append(unmapped_spec)

    return formatted, is_qcfail


# Do the mappy alignment
def mappy_al(rcon: str, rfmt: str, cpus: int, tech: str, ref: str, outf: str, foutf:str, file1:list, file2:list=None):
    """Align single-end reads to reference.
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
#    aligner = mappy_al_ref(ref, tech, cpus) # create aligner
    flagstats = {'total_reads':0, 'unmapped':0,'removed_reads_primary':0, 'kept_primary':0, 'kept_secondary':0, 'kept_supplementary':0} # create flagstats dictionary
    if rcon == "paired-end":
        flagstats['read_1'] = 0
        flagstats['read_2'] = 0
    sq_head, read_count, read_iter = alignment_preprocessing(ref, rfmt, file1) if rcon == "single-end" else alignment_preprocessing(ref, rfmt, file1, file2) # create header and get totalcount of reads
    obam = pysam.AlignmentFile(outf, 'wb', header=sq_head) # create output BAM file
    fbam = pysam.AlignmentFile(foutf, 'wb', header=sq_head) # create failed BAM file
    flagstats['total_reads'] = read_count
    # for scampiman progress bar
    pool_tally = 0
    chunk_size = 10000 if int(read_count / cpus) < 10000 else int(read_count / cpus) # after how many reads to update progress bar
    progress_list = [read_num * chunk_size for read_num in range(0,int(read_count/chunk_size))] # create list of read counts to update progress bar based on multiplier
    progress_list.append(read_count) # add total number of reads to progress list
    header_starttime = time.perf_counter() # start timer for progress bar
    
    worker_fn = align_se_chunk if rcon == "single-end" else align_pe_chunk

    # Dispatch chunks to the process pool and write results as they arrive
    with ProcessPoolExecutor(
        max_workers=cpus,
        initializer=mappy_al_ref,
        initargs=(ref, tech),
    ) as executor:

        futures = {
            executor.submit(worker_fn, chunk): idx
            for idx, chunk in enumerate(_chunked(read_iter, chunk_size))
        }

        for future in as_completed(futures):
            # Header stuff
            header_endtime = time.perf_counter()
            time_taken = header_endtime - header_starttime
            shrimp_progress(read_count, progress_list[pool_tally], time_taken, "align")
            pool_tally += 1

            for result in future.result():
                if not result['hits']:
                    rec = pysam.AlignedSegment(header=sq_head)
                    rec.query_name      = result['query_name']
                    rec.query_sequence  = result['query_sequence']
                    rec.query_qualities = result['query_qualities']
                    rec.is_unmapped     = True
                    flagstats['unmapped'] += 1

                    if rfmt == "bam":
                        rec.tags = result['original_tags']
                    if rcon == "paired-end":
                        rec2 = pysam.AlignedSegment(header=sq_head)
                        rec2.query_name      = result['query_name2']
                        rec2.query_sequence  = result['query_sequence2']
                        rec2.query_qualities = result['query_qualities2']
                        rec2.is_unmapped     = True
                        flagstats['unmapped'] += 1
                        fbam.write(rec2)

                    fbam.write(rec)
                    continue

                if result['is_qcfail']:
                    flagstats['removed_reads_primary'] += 1 if rcon == "single-end" else 2

                if rcon == "single-end":
                    original_tags = result['original_tags']
                    sa_tags = [fmt['sa_tag_str'] for fmt in result['formatted']]

                # write to output BAM
                for hit_index, fmt in enumerate(result['formatted']):
                    rec = pysam.AlignedSegment(header=sq_head)
                    rec.query_name        = fmt['query_name']
                    rec.query_sequence    = fmt['seq']
                    rec.query_qualities   = fmt['qual']
                    rec.reference_name    = fmt['reference_name']
                    rec.reference_start   = fmt['reference_start']
                    rec.mapping_quality   = fmt['mapping_quality']
                    rec.cigarstring       = fmt['cigarstring']
                    rec.is_forward        = fmt['is_forward']
                    rec.is_reverse        = fmt['is_reverse']
                    rec.is_supplementary  = fmt['is_supplementary']
                    rec.is_secondary      = fmt['is_secondary']
                                            
                    if rcon == "paired-end":
                        rec.is_read1      = fmt['is_read1']
                        rec.is_read2      = fmt['is_read2']
                        rec.is_paired         = fmt['is_paired']
                        rec.is_proper_pair    = fmt['is_proper_pair']
                        rec.next_reference_name = fmt['next_reference_name']
                        rec.next_reference_start = fmt['next_reference_start']
                        rec.template_length = fmt['template_length']
                        rec.mate_is_forward   = fmt['mate_is_forward']
                        rec.mate_is_reverse   = fmt['mate_is_reverse']
                        rec.tags = fmt['tags']
                    else: # single-end
                        read_sa_tag = sa_tags.copy()
                        read_sa_tag.pop(hit_index)
                        base_tags = original_tags + fmt['tags'] if rfmt == "bam" else fmt['tags']
                        rec.tags = base_tags + [("SA", ';'.join(read_sa_tag))] if read_sa_tag else base_tags

                    if result['is_qcfail']:
                        fbam.write(rec)
                        continue

                    # Add in flagstats counting
                    if fmt['is_secondary']:
                        flagstats['kept_secondary'] += 1    
                    elif fmt['is_supplementary']:
                        flagstats['kept_supplementary'] += 1
                    else:
                        flagstats['kept_primary'] += 1
                    if rcon == "paired-end":
                        if fmt['is_read1']:
                            flagstats['read_1'] += 1
                        elif fmt['is_read2']:
                            flagstats['read_2'] += 1

                    obam.write(rec)

    obam.close()
    fbam.close() 
    pysam.sort("-o", outf.replace('.bam', '.sort.bam'), outf) # sort BAM file
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

