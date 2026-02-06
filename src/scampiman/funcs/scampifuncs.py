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


def mappy_al_ref(ref: str, tech: str, cpus:int):
    if tech == "ont":
        targ = "lr:hq"
    elif tech == "illumina":
        targ = "sr"

    aligner = mp.Aligner(ref,
        preset=targ, n_threads = cpus)
    
    return aligner


def mappy_al_header(ref: str, rfmt:str, file:str):
    if rfmt == "bam":
        sq_head = pysam.AlignmentFile(file, 'rb', check_sq=False).header.to_dict() | \
            pysam.AlignmentHeader(
            ).from_references(
                reference_names=pysam.FastaFile(ref).references, 
                reference_lengths=pysam.FastaFile(ref).lengths
                ).to_dict()
    elif rfmt == "fastq":
        sq_head = pysam.AlignmentHeader(
            ).from_references(
                reference_names=pysam.FastaFile(ref).references, 
                reference_lengths=pysam.FastaFile(ref).lengths
                ).to_dict()
    return sq_head



def mappy_hits_bam_fmt_single(q: dict, hits: list, header_dict: dict):
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
        
        if str(hit) == str(hits[0]):
            if hit.strand == 1: 
                rec.is_forward = True
                rec.query_sequence = q['query_sequence']
                rec.query_qualities = q['query_qualities']
                cigar_start = str(hit.q_st)
                cigar_end = str(len(q['query_sequence'])-hit.q_en)
            elif hit.strand == -1: 
                rec.is_reverse = True
                rec.query_sequence = mp.revcomp(q['query_sequence'])
                rec.query_qualities = q['query_qualities'][::-1]
                cigar_start = str(len(q['query_sequence'])-hit.q_en)
                cigar_end = str(hit.q_st)
            
            if cigar_start != '0' and cigar_end != '0':
                rec.cigarstring = cigar_start + 'S' + hit.cigar_str + cigar_end + 'S'
            elif cigar_start != '0':
                rec.cigarstring = cigar_start + 'S' + hit.cigar_str
            elif cigar_end != '0':
                rec.cigarstring = hit.cigar_str + cigar_end + 'S'
            else:
                rec.cigarstring = hit.cigar_str 
            recs_list.append(rec)
        else:
            prime_qrange = list(range(hits[0].q_st,hits[0].q_en))
            hit_qrange = list(range(hit.q_st, hit.q_en))
            common_bases = set(prime_qrange) & set(hit_qrange)
            qlen_olap_ratio = len(common_bases)/len(hit_qrange)
            if qlen_olap_ratio < 0.5: 
                rec.is_supplementary = True
            else:
                rec.is_secondary = True

            if hit.strand == 1: 
                rec.is_forward = True
                rec.query_sequence = q['query_sequence'][hit.q_st:hit.q_en]
                rec.query_qualities = q['query_qualities'][hit.q_st:hit.q_en]
                rec.cigarstring = str(hit.q_st) + 'H' + hit.cigar_str + str(len(q['query_sequence'])-hit.q_en) + 'H'
            if hit.strand == -1: 
                rec.is_reverse = True
                rec.query_sequence = mp.revcomp(q['query_sequence'][hit.q_st:hit.q_en])
                rec.query_qualities = q['query_qualities'][hit.q_st:hit.q_en][::-1]
                rec.cigarstring =  str(len(q['query_sequence'])-hit.q_en) +'H' + hit.cigar_str + str(hit.q_st) + 'H'

            ref_olap = recs_list[0].get_overlap(hit.r_st, hit.r_en) 
            hit_refRange = hit.r_en - hit.r_st
            rlen_olap_ratio = ref_olap/hit_refRange
            ## if the overlap is less than 50% of the hit, it is a secondary alignment
            if rlen_olap_ratio < 0.5:
                #rec.is_qcfail = True
                recs_list[0].is_qcfail = True
            elif rec.is_supplementary and hits[0].strand == hit.strand:
                #rec.is_qcfail = True
                recs_list[0].is_qcfail = True

            recs_list.append(rec)
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

        if not hit.is_primary:
            tp = 'S'
        else:
            tp = 'P'

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
        rec.query_name = q['query_name']
        rec.reference_name = hit.ctg 
        rec.reference_start = hit.r_st 
        rec.mapping_quality = hit.mapq
        rec.is_paired = True
        if not hit.is_primary:
            rec.is_secondary = True
            tp = 'S'
        else:
            tp = 'P'
        
        if hit.read_num == 1:
            rec.is_read1 = True
            q_seq = q['seq1']
            q_qual = q['quality1']
            hit_r1.append(hit)
        elif hit.read_num == 2:
            rec.is_read2 = True
            q_seq = q['seq2']
            q_qual = q['quality2']
            hit_r2.append(hit)

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

    return rec_list




def mappy_al(rfmt: str, cpus: int, tech: str, ref: str,  outf: str, foutf:str, file1:str, file2:str=None):
    aligner = mappy_al_ref(ref, tech, cpus)
    flagstats = {'total_reads':0, 'unmapped':0,'removed_reads_primary':0, 'kept_primary':0, 'kept_secondary':0, 'kept_supplementary':0}
    sq_head = mappy_al_header(ref, rfmt, file1)
    obam = pysam.AlignmentFile(outf, 'wb', header=sq_head)
    fbam = pysam.AlignmentFile(foutf, 'wb', header=sq_head)
    
    if rfmt == "bam":
        if tech == "ont":
            for i in pysam.AlignmentFile(file1, 'rb', check_sq=False):
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
        if file2 and tech == "illumina":
            flagstats = flagstats | {'read1':0, 'read2':0}
            with pysam.FastxFile(file1) as fq1, pysam.FastxFile(file2) as fq2:
                for read1, read2 in zip(fq1, fq2):
                    flagstats['total_reads'] += 2
                    hits = list(aligner.map(read1.sequence, seq2=read2.sequence, cs=True, MD=True))
                    read_name = list(set([read1.name.split("/")[0], read2.name.split("/")[0]]))[0]

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

                    rec_list = mappy_hits_bam_fmt_paired({"query_name": read_name, "seq1": read1.sequence, "seq2": read2.sequence, "quality1": read1.quality, "quality2": read2.quality}, hits, sq_head)

                    for rec in rec_list:
                        flagstats['kept_primary'] += 1
                        if rec.is_read1:
                            flagstats['read1'] += 1
                        elif rec.is_read2:
                            flagstats['read2'] += 1
                        obam.write(rec)

        else:
            for i in pysam.FastxFile(file1):
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

    pysam.sort("-o", outf.replace('.bam', '.sorted.bam'), outf)
    pysam.index("-b", "-o", outf.replace('.bam', '.sorted.bam') + '.bai', outf.replace('.bam', '.sorted.bam'))
   
    return flagstats



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


def flag_stats(bam: str, outf: str):

    flagstats = open(outf, 'a')
    flagstats_command = ['samtools', 'flagstats', '-O', 'tsv', bam]
    logger.info(flagstats_command)

    # Use subprocess to capture the STDOUT
    flagstats_process = subprocess.Popen(
        flagstats_command, 
        stdout=flagstats
    )

    output, _ = flagstats_process.communicate()
    if output:
        return output
    else:
        return None
