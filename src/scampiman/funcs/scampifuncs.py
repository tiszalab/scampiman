import logging
import re, gzip, shutil
import argparse
import pysam
import sys, os, time
import pandas as pd
import mappy as mp
from typing import Optional, Tuple, List, Dict, Any
from datetime import timedelta

# since I'm logging, I need the same logger as in the main script
logger = logging.getLogger("pct_logger")

def pysam_cat(reads: str, ftype: str, outf: str):

    if ftype == "directory":
        bam_list = []
        for directory in reads.split():
            for bam in os.listdir(directory):
                if bam.endswith('.bam'):
                    f = os.path.join(directory, bam)

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

def shrimp_progress(total_reads: int, aligned_reads: int, time_taken: float, bar_width: int) -> None:
    """Display a progress bar with a shrimp swimming across the terminal.
    
    Args:
        total_reads: Total number of reads to align.
        bar_width: Width of the progress bar in characters.
    """
    if aligned_reads == 0:
        print(f"\n🌊 Aligning {total_reads} reads with shrimp power! 🌊\n")

    progress = aligned_reads / total_reads
    shrimp_pos = int(progress * bar_width)
    
    # Build the progress bar
    water_before = "~" * shrimp_pos
    water_after = "~" * (bar_width - shrimp_pos)
    
    # Create the bar with shrimp at current position
    bar = f"[{water_before}🦐{water_after}]"
    
    # Calculate percentage
    percent = int(progress * 100)
    # Print the progress bar (overwrite previous line)
    sys.stdout.write(f"\r{bar} {percent:3d}% | {aligned_reads}/{total_reads} | {timedelta(seconds=time_taken)}")
    sys.stdout.flush()
    
    if aligned_reads == total_reads - 1:
        print("\n\n✨ Shrimp has arrived! Alignment complete! ✨\n")
    elif aligned_reads < total_reads:
        time.sleep(1)



def file_paths(reads: str, rfmt: str, rcon: str, intype: str):
    """ Makes a list of paths to  input files from a string. 
    This is mostly for finding the files within a directory and/or return matching R1 and R2 files.
    If intype is "files" and the rcon is "single-end", it will just return a list version of the "reads" string.

    Args:
        reads: string of space separated paths to input files or directories
        rfmt: "bam" or "fastq"
        rcon: "single-end" or "paired-end"
        intype: "directory" or "files"
    """
    file_list = []
    # establish the file extensions
    if rfmt == "bam":
        extensions = ('.bam')
    elif rfmt == "fastq":
        extensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
    # find the files with file extensions in directories
    if intype == "directory":
        for directory in reads.split():
            for file in os.listdir(directory):
                if file.endswith(extensions):
                    f = os.path.join(directory, file)
                    if os.path.isfile(f) and os.path.getsize(f) > 0:
                        file_list.append(f)
    # if string of files given, filter for file extensions and put into list
    elif intype == "files": 
        for file in reads.split():
            if file.endswith(extensions):
                f = os.path.abspath(file)
                if os.path.isfile(f) and os.path.getsize(f) > 0:
                    file_list.append(f)
                else:
                    logger.info(f"File {file} excluded; Does not exist or is empty.") 
            else:
                logger.info(f"File {file} excluded; Does not end with {extensions}") 
    # if single-end, return the list of files
    if rcon == "single-end":
        logger.info(f"Single-end files: {file_list}")
        return file_list 
    # if paired-end, sort the files into R1 and R2
    elif rcon == "paired-end":
        #establish the file patterns
        r1_pattern = re.compile(r'(.+?)(_R1_|_R1\.|.R1\.|_1\.|_1_)(.*)$')
        r2_pattern = re.compile(r'(.+?)(_R2_|_R2\.|.R2\.|_2\.|_2_)(.*)$')
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
        for key in r1_files.keys():
            if key in r2_files.keys():
                matched_r1.append(r1_files[key])
                matched_r2.append(r2_files[key])
            else:
                logger.info(f"No R2 pair found for: {r1_files[key]}") # if it doesn't have a matching R2 file, print a message
        # look for R2 files that don't have a matching R1 file
        for key in r2_files.keys():
            if key not in r1_files.keys():
                logger.info(f"No R1 pair found for: {r2_files[key]}") # print message  
        # return the matched read1 and read2 lists
        logger.info(f"Paired-end files: {matched_r1, matched_r2}")
        return matched_r1, matched_r2



def mappy_al_ref(ref: str, tech: str, cpus:int):
    if tech == "ont":
        targ = "lr:hq"
    elif tech == "illumina":
        targ = "sr"

    aligner = mp.Aligner(ref,
        preset=targ, n_threads = cpus)
    
    return aligner


def mappy_al_header(ref: str, rfmt:str, file_list: list):
    """Make a header for the alignment file."""
    read_count = 0
    sq_head = pysam.AlignmentHeader(
    ).from_references(
        reference_names=pysam.FastaFile(ref).references, 
        reference_lengths=pysam.FastaFile(ref).lengths
        ).to_dict()
    sq_head['PG'] = []
    if rfmt == "bam":
        for file in file_list:
            #bam = pysam.AlignmentFile(file, 'rb', check_sq=False)
            read_count += pysam.AlignmentFile(file, 'rb', check_sq=False).count(until_eof=True)
            if len(sq_head) == 2: # the only thing in the freshly made header is the SQ tag and an empty PG tag
                sq_head = pysam.AlignmentFile(file, 'rb', check_sq=False).header.to_dict() | sq_head # combine the headers
            else: # if the header has already been combined
                sq_head['RG'] += pysam.AlignmentFile(file, 'rb', check_sq=False).header.to_dict()['RG'] # add in the RG tag just incase smaples were from different runs
        # get rid of duplicate RG tags from looping
        g = list({frozenset(d.items()) for d in sq_head['RG']}) # convert the list into a set to remove duplicates
        g = [dict(f) for f in g] # Convert back to dictionaries
        sq_head['RG'] = g # set RG header tag
    elif rfmt == "fastq":
        for file in file_list:
            for i in pysam.FastxFile(file):
                read_count += 1

    # add in the PG tag (gives details about mappy and which read files were used for alignment)
    sq_head['PG'].append({'ID': 'aligner', 'PN': 'mappy', 'VN': mp.__version__, 'DS': 'minimap2 alignment', 'CL': ' '.join(file_list)})

    return sq_head, read_count



def mappy_hits_bam_fmt_single(q: dict, hits: list, header_dict: dict):
    """Format mappy hits to bam format
    
    Args:
        q: query dictionary (read name, sequence, quality/qscore string)
        hits: list of mappy hits
        header_dict: header dictionary from mappy_al_header or pysam.AlignmentHeader.from_dict()
    """
    sq_head = pysam.AlignmentHeader.from_dict(header_dict)
    recs_list = []
    sa_tags = []
    tags = []

    for hit in hits:
        rec = pysam.AlignedSegment(header=sq_head)
        rec.is_mapped = True
        rec.query_name = q['query_name']
        rec.reference_name = hit.ctg 
        rec.reference_start = hit.r_st 
        rec.mapping_quality = hit.mapq 
        
        if str(hit) == str(hits[0]): # only set these tags for the primary hit (str(hits[0]) is the primary hit)
            if hit.strand == 1: # forward strand
                rec.is_forward = True
                rec.query_sequence = q['query_sequence']
                rec.query_qualities = q['query_qualities']
                cigar_front_clip = str(hit.q_st) # which base of the query sequence is the first base of the alignment; how many bases are clipped from the front
                cigar_end_clip = str(len(q['query_sequence'])-hit.q_en) # which base of the query sequence is the last base of the alignment; how many bases are clipped from the end
            elif hit.strand == -1: # reverse strand
                rec.is_reverse = True
                rec.query_sequence = mp.revcomp(q['query_sequence']) # bam output for reverse strand is reverse complement
                rec.query_qualities = q['query_qualities'][::-1] # the quality string is also reverse
                cigar_front_clip = str(len(q['query_sequence'])-hit.q_en) # which base of the query sequence (as reverse) is the first base of the alignment
                cigar_end_clip = str(hit.q_st) # which base of the query sequence (as reverse) is the last base of the alignment
            # build cigar string
            if cigar_front_clip != '0' and cigar_end_clip != '0': # if there are clippings at both ends
                rec.cigarstring = cigar_front_clip + 'S' + hit.cigar_str + cigar_end_clip + 'S'
            elif cigar_front_clip != '0': # if there is a clipping at the front but not the back
                rec.cigarstring = cigar_front_clip + 'S' + hit.cigar_str
            elif cigar_end_clip != '0': # if there is a clipping at the end but not the front
                rec.cigarstring = hit.cigar_str + cigar_end_clip + 'S'
            else: # if there are no clippings
                rec.cigarstring = hit.cigar_str 
            recs_list.append(rec)
        else:
            prime_qrange = list(range(hits[0].q_st,hits[0].q_en)) # range of base positions on the read of the primary read alignment hits[0]
            hit_qrange = list(range(hit.q_st, hit.q_en)) # range of base positions on the read of the hit
            common_bases = set(prime_qrange) & set(hit_qrange) # common READ bases between the primary hit and the hit
            qlen_olap_ratio = len(common_bases)/len(hit_qrange) # ratio of common READ bases to the length of the hit
            if qlen_olap_ratio < 0.5: 
                rec.is_supplementary = True
            else:
                rec.is_secondary = True

            if hit.strand == 1: 
                rec.is_forward = True
                rec.query_sequence = q['query_sequence'][hit.q_st:hit.q_en] # report only the aligned sequence bases
                rec.query_qualities = q['query_qualities'][hit.q_st:hit.q_en] # report only the aligned quality positions
                rec.cigarstring = str(hit.q_st) + 'H' + hit.cigar_str + str(len(q['query_sequence'])-hit.q_en) + 'H' # building cigar string for secondary alignment; H = hard clipping
            if hit.strand == -1: 
                rec.is_reverse = True
                rec.query_sequence = mp.revcomp(q['query_sequence'][hit.q_st:hit.q_en]) # report only the aligned sequence bases as reverse complement
                rec.query_qualities = q['query_qualities'][hit.q_st:hit.q_en][::-1] # report only the aligned quality positions as reverse
                rec.cigarstring =  str(len(q['query_sequence'])-hit.q_en) +'H' + hit.cigar_str + str(hit.q_st) + 'H' # building cigar string for secondary alignment; H = hard clipping

            ref_olap = recs_list[0].get_overlap(hit.r_st, hit.r_en) # number of overlapping reference bases between the primary hit and the hit
            hit_refRange = hit.r_en - hit.r_st # reference length of the hit
            rlen_olap_ratio = ref_olap/hit_refRange # ratio of overlapping reference bases of the primary and the hit to the length of the hit; how much of the hit is covered by the primary hit
            if rlen_olap_ratio < 0.5: ## if the overlap is less than 50% of the hit, it is a secondary alignment 
                recs_list[0].is_qcfail = True # set primary alignment as qcfail to filter out read later
            elif rec.is_supplementary and hits[0].strand == hit.strand: # if the hit is supplementary and on the same strand as the primary hit
                recs_list[0].is_qcfail = True # set primary alignment as qcfail to filter out read later

            recs_list.append(rec)
        
        # if the primary alignment is qcfail, skip the tag creations
        if recs_list[0].is_qcfail:
            continue
        ## creating tags
        ### de tag: Gap-compressed per-base sequence divergence = matching 
        de_denom = 0 #the de denominator
        for i in hit.cigar: # loop through the array of CIGAR [length, operation]
            if i[1] == 0: #operation 0: Number of matched bases; Add up the number of matched bases
                de_denom += i[0] # Add up the number of matched bases for the denominator
            else:
                de_denom += 1 # Add 1 for any other operation (insertion, deletion, etc.)
        de = 1-(hit.mlen/de_denom)
        rounded_de = round(de, str(round(de,6)).split('.')[1].count('0') + 6)
        ### tp tag: primary (P) or [supplementary or secondary] (S)
        if not hit.is_primary:
            tp = 'S'
        else:
            tp = 'P'
        ### strand char for SA tag  
        if hit.strand == 1: 
            strand_char = '+'
        elif hit.strand == -1:
            strand_char = '-'

        sa_tags.append(f"{hit.ctg},{hit.r_st},{strand_char},{rec.cigarstring},{hit.mapq},{hit.NM}")
        tags.append([("NM", hit.NM), ("nn", (hit.blen - hit.mlen - hit.NM), "i"),("tp", tp, "A"),("de", rounded_de, "f"),("MD", hit.MD), ("cs", hit.cs)])
    
    return recs_list, tags, sa_tags


def mappy_hits_bam_fmt_paired(q: dict, hits: list, header_dict: dict):
    sq_head = pysam.AlignmentHeader.from_dict(header_dict)
    rec_list = []
    hit_r1 = []
    hit_r2 = []
    for hit in hits:
        rec = pysam.AlignedSegment(header=sq_head)
        rec.is_mapped = True
        rec.reference_name = hit.ctg 
        rec.reference_start = hit.r_st 
        rec.mapping_quality = hit.mapq
        rec.is_paired = True
        
        if hit.read_num == 1:
            rec.is_read1 = True
            rec.query_name = f"{q['query_name_read1']}_R1"
            q_seq = q['seq1']
            q_qual = q['quality1']
            hit_r1.append(hit)

        elif hit.read_num == 2:
            rec.is_read2 = True
            rec.query_name = f"{q['query_name_read2']}_R2"
            q_seq = q['seq2']
            q_qual = q['quality2']
            hit_r2.append(hit)

        if not hit.is_primary:
            tp = 'S'
            if hit.read_num == 1:
                prime_qrange = list(range(hit_r1[0].q_st,hit_r1[0].q_en))
            elif hit.read_num == 2:
                prime_qrange = list(range(hit_r2[0].q_st,hit_r2[0].q_en))
            
            hit_qrange = list(range(hit.q_st, hit.q_en))
            common_bases = set(prime_qrange) & set(hit_qrange)
            qlen_olap_ratio = len(common_bases)/len(hit_qrange)
            
            if qlen_olap_ratio < 0.5: 
                rec.is_supplementary = True
            else:
                rec.is_secondary = True     
        else:
            tp = 'P'

        if hit.strand == 1:
            rec.is_forward = True
            rec.mate_is_reverse = True
            rec.query_sequence = q_seq
            rec.query_qualities = q_qual
            cigar_start = str(hit.q_st)
            cigar_end = str(len(q_qual)-hit.q_en)
        elif hit.strand == -1: 
            rec.is_reverse = True
            rec.mate_is_forward = True
            rec.query_sequence = mp.revcomp(q_seq)
            rec.query_qualities = q_qual[::-1]
            cigar_start = str(len(q_qual)-hit.q_en)
            cigar_end = str(hit.q_st)
        
        if cigar_start != '0' and cigar_end != '0':
            rec.cigarstring = cigar_start + 'S' + hit.cigar_str + cigar_end + 'S'
        elif cigar_start != '0':
            rec.cigarstring = cigar_start + 'S' + hit.cigar_str
        elif cigar_end != '0':
            rec.cigarstring = hit.cigar_str + cigar_end + 'S'
        else:
            rec.cigarstring = hit.cigar_str
        ## creating tags
        ### de tag: Gap-compressed per-base sequence divergence = matching 
        de_denom = 0 #the de denominator
        for i in hit.cigar: # loop through the array of CIGAR [length, operation]
            #cig_dict[i[1]] += i[0] #for every operation, add the length to the dictionary
            if i[1] == 0: #operation 0: Number of matched bases; Add up the number of matched bases
                de_denom += i[0] # Add up the number of matched bases for the denominator
            else:
                de_denom += 1 # Add 1 for any other operation (insertion, deletion, etc.)
        de = 1-(hit.mlen/de_denom)
        rounded_de = round(de, str(round(de,6)).split('.')[1].count('0') + 6)
        rec.tags =  [("NM", hit.NM), ("nn", (hit.blen - hit.mlen - hit.NM), "i"),("tp", tp, "A"), ("de", rounded_de, "f"),("MD", hit.MD), ("cs", hit.cs)]

        rec_list.append(rec)

    for h1, h2 in zip(hit_r1, hit_r2):
        rec_1 = rec_list[hits.index(h1)]
        rec_2 = rec_list[hits.index(h2)]
        rec_1.next_reference_name = rec_2.reference_name
        rec_1.next_reference_id = rec_2.reference_id
        rec_1.next_reference_start = rec_2.reference_start
        rec_1.is_proper_pair = True
        rec_2.next_reference_id = rec_1.reference_id
        rec_2.next_reference_start = rec_1.reference_start
        rec_2.next_reference_name = rec_1.reference_name
        rec_2.is_proper_pair = True
        
        if not rec_1.is_reverse:
            rec_2.template_length = h2.r_st - h1.r_st + h2.mlen
            rec_1.template_length = rec_2.template_length * -1
        else:
            rec_1.template_length = h1.r_st - h2.r_st + h1.mlen
            rec_2.template_length = rec_1.template_length * -1
        # add them back in
        rec_list[hits.index(h1)] = rec_1
        rec_list[hits.index(h2)] = rec_2

        if rec_1.is_forward == rec_2.is_forward: #if both reads are on the same strand
            rec_list[0].is_qcfail = True
        if rec_1.get_overlap(h2.r_st, h2.r_en) == 0: # if the reference mapping is not overlapping between the two reads
            rec_list[0].is_qcfail = True
    
    if len(hit_r1) != len(hit_r2): #if there are not the same number of hits for both reads
        rec_list[0].is_qcfail = True
    # add in unmapped reads to rec_list to be add into the fail BAM file
    if not hit_r1 or not hit_r2:
        rec_list[0].is_qcfail = True
        unmapped_read = pysam.AlignedSegment(header=sq_head)
        unmapped_read.is_unmapped = True
        unmapped_read.is_paired = True
        if not hit_r1:
            unmapped_read.is_read1 = True
            unmapped_read.query_name = f"{q['query_name_read1']}_R1"
            unmapped_read.query_sequence = q['seq1']
            unmapped_read.query_qualities = q['quality1']
        elif not hit_r2:
            unmapped_read.is_read2 = True
            unmapped_read.query_name = f"{q['query_name_read2']}_R2"
            unmapped_read.query_sequence = q['seq2']
            unmapped_read.query_qualities = q['quality2']
        rec_list.append(unmapped_read)

    return rec_list



def mappy_al_single(rfmt: str, cpus: int, tech: str, ref: str, outf: str, foutf:str, file1:list):
    aligner = mappy_al_ref(ref, tech, cpus)
    flagstats = {'total_reads':0, 'unmapped':0,'removed_reads_primary':0, 'kept_primary':0, 'kept_secondary':0, 'kept_supplementary':0}
    sq_head, read_count = mappy_al_header(ref, rfmt, file1)
    obam = pysam.AlignmentFile(outf, 'wb', header=sq_head)
    fbam = pysam.AlignmentFile(foutf, 'wb', header=sq_head)
    
    # for scampiman progress bar
    multiplier = int(read_count / 50)
    progress_list = [read_count * multiplier for read_count in range(0,50)]
    progress_list.append(read_count - 1)
    header_starttime = time.perf_counter()
    if rfmt == "bam":
        for bam in file1:
            for i in pysam.AlignmentFile(bam, 'rb', check_sq=False):
                if flagstats['total_reads'] in progress_list:
                    header_endtime = time.perf_counter()
                    time_taken = header_endtime - header_starttime

                    
                    shrimp_progress(read_count, flagstats['total_reads'], time_taken, 40)
            
                flagstats['total_reads'] += 1
                hits = list(aligner.map(i.query_sequence, cs=True, MD=True))
                if not hits:
                    i.is_unmapped = True
                    flagstats['unmapped'] += 1
                    fbam.write(i)
                    continue

                recs_list, tags, sa_tags = mappy_hits_bam_fmt_single({"query_name": i.query_name, "query_sequence": i.query_sequence, "query_qualities": i.query_qualities}, hits, sq_head)
            
                if recs_list[0].is_qcfail:
                    flagstats['removed_reads_primary'] += 1
                    for rec in recs_list:
                        fbam.write(rec)
                    continue        
                
                for read_alignment in recs_list:
                    if read_alignment.is_secondary:
                        flagstats['kept_secondary'] += 1
                    elif read_alignment.is_supplementary:
                        flagstats['kept_supplementary'] += 1
                    else:
                        flagstats['kept_primary'] += 1

                    hit_index = recs_list.index(read_alignment)
                    read_sa_tag = sa_tags.copy()
                    read_sa_tag.pop(hit_index)
                    if not read_sa_tag:
                        read_alignment.tags =  i.get_tags() + tags[hit_index]
                        obam.write(read_alignment)
                        continue
                    
                    read_alignment.tags =  i.get_tags() + tags[hit_index] + [("SA", ';'.join(read_sa_tag))]
                    obam.write(read_alignment)

    elif rfmt == "fastq":
        for fastq in file1:
            for i in pysam.FastxFile(fastq):
                if flagstats['total_reads'] in progress_list:
                    shrimp_progress(read_count, flagstats['total_reads'], 40)
                flagstats['total_reads'] += 1
                hits = list(aligner.map(i.sequence, cs=True, MD=True))
                if not hits:
                    rec = pysam.AlignedSegment(header=pysam.AlignmentHeader.from_dict(sq_head))
                    rec.query_sequence = i.sequence
                    rec.query_qualities = i.quality
                    rec.query_name = i.name
                    rec.is_unmapped = True
                    flagstats['unmapped'] += 1
                    fbam.write(rec)
                    continue

                recs_list, tags, sa_tags = mappy_hits_bam_fmt_single({"query_name": i.name, "query_sequence": i.sequence, "query_qualities": i.quality}, hits, sq_head)
            
                if recs_list[0].is_qcfail:
                    flagstats['removed_reads_primary'] += 1
                    for rec in recs_list:
                        fbam.write(rec)
                    continue        
                
                for read_alignment in recs_list:
                    if read_alignment.is_secondary:
                        flagstats['kept_secondary'] += 1
                    elif read_alignment.is_supplementary:
                        flagstats['kept_supplementary'] += 1
                    else:
                        flagstats['kept_primary'] += 1

                    hit_index = recs_list.index(read_alignment)
                    read_sa_tag = sa_tags.copy()
                    read_sa_tag.pop(hit_index)
                    if not read_sa_tag:
                        read_alignment.tags = tags[hit_index]
                        obam.write(read_alignment)
                        continue
                    
                    read_alignment.tags = tags[hit_index] + [("SA", ';'.join(read_sa_tag))]
                    obam.write(read_alignment)
    obam.close()
    fbam.close() 
    pysam.sort("-o", outf.replace('.bam', '.sort.bam'), outf)
    
    flag_df = pd.DataFrame(flagstats, index=[0])
    return flag_df



def mappy_al_paired(cpus: int, tech: str, ref: str, outf: str, foutf:str, file1:list, file2:list):
    aligner = mappy_al_ref(ref, tech, cpus)
    flagstats = {'total_reads':0, 'unmapped':0,'removed_reads_primary':0, 'kept_primary':0, 'kept_secondary':0, 'kept_supplementary':0, 'read1':0, 'read2':0}

    files = file1 + file2
    sq_head, read_count = mappy_al_header(ref, rfmt, files)

    # for scampiman progress bar
    multiplier = int(read_count / 40)
    progress_list = [read_count * multiplier for read_count in range(0,40)]
    progress_list.append(read_count)

    obam = pysam.AlignmentFile(outf, 'wb', header=sq_head)
    fbam = pysam.AlignmentFile(foutf, 'wb', header=sq_head)

    for fastq1, fastq2 in zip(file1, file2):
        with pysam.FastxFile(fastq1) as fq1, pysam.FastxFile(fastq2) as fq2:
            for read1, read2 in zip(fq1, fq2):
                if flagstats['total_reads'] in progress_list:
                    shrimp_progress(read_count, flagstats['total_reads'], 40)
                flagstats['total_reads'] += 2

                hits = list(aligner.map(read1.sequence, seq2=read2.sequence, cs=True, MD=True))
                if not hits:
                    rec1 = pysam.AlignedSegment(header=pysam.AlignmentHeader.from_dict(sq_head))
                    rec2 = pysam.AlignedSegment(header=pysam.AlignmentHeader.from_dict(sq_head))
                    rec1.query_sequence = read1.sequence
                    rec2.query_sequence = read2.sequence
                    rec1.query_qualities = read1.quality
                    rec2.query_qualities = read2.quality
                    rec1.is_read1 = True
                    rec2.is_read2 = True
                    for i in [rec1, rec2]:
                        i.query_name = read_name
                        i.is_unmapped = True
                        i.is_paired = True
                    fbam.write(rec1)
                    fbam.write(rec2)
                    flagstats['unmapped'] += 2
                    continue

                rec_list = mappy_hits_bam_fmt_paired({"query_name_read1": read1.name, "query_name_read2": read2.name,"seq1": read1.sequence, "seq2": read2.sequence, "quality1": read1.quality, "quality2": read2.quality}, hits, sq_head)

                if rec_list[0].is_qcfail:
                    flagstats['removed_reads_primary'] += 2
                    for rec in rec_list:
                        print(rec)
                        fbam.write(rec)
                    continue
                    
                for rec in rec_list:                            
                    flagstats['kept_primary'] += 1
                    if rec.is_read1:
                        flagstats['read1'] += 1
                    elif rec.is_read2:
                        flagstats['read2'] += 1
                    obam.write(rec)               
    obam.close()
    fbam.close() 
    pysam.sort("-o", outf.replace('.bam', '.sort.bam'), outf)
    
    flag_df = pd.DataFrame(flagstats, index=[0])
    return flag_df



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
    #max_length += 50
    return max_length * 2

