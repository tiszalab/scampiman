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


def mappy_al_ref(ref: str, tech: str, cpus:int):
    if tech == "ont":
        targ = "lr:hq"
    elif tech == "illumina":
        targ = "sr"

    aligner = mp.Aligner(ref,
        preset=targ, n_threads = cpus)
    
    return aligner


def mappy_hits_bam_fmt(q: pysam.libcalignedsegment.AlignedSegment, hits: list, header_dict: dict):
    sq_head = pysam.AlignmentHeader.from_dict(header_dict)
    recs_list = []
    sa_tags = []
    tags = []

    for hit in hits:
        rec = pysam.AlignedSegment(header=sq_head)
        rec.is_mapped = True
        rec.query_name = q.query_name
        rec.reference_name = hit.ctg 
        rec.reference_start = hit.r_st 
        rec.mapping_quality = hit.mapq 
        
        if str(hit) == str(hits[0]):
            if hit.strand == 1: 
                rec.is_forward = True
                rec.query_sequence = q.query_sequence
                rec.query_qualities = q.query_qualities
                rec.cigarstring = str(hit.q_st) + 'S' + hit.cigar_str + str(len(q.query_sequence)-hit.q_en) + 'S'
            if hit.strand == -1: 
                rec.is_reverse = True
                rec.query_sequence = mp.revcomp(q.query_sequence)
                rec.query_qualities = q.query_qualities[::-1]
                rec.cigarstring =  str(len(q.query_sequence)-hit.q_en) +'S' + hit.cigar_str + str(hit.q_st) + 'S'
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
                rec.query_sequence = q.query_sequence[hit.q_st:hit.q_en]
                rec.query_qualities = q.query_qualities[hit.q_st:hit.q_en]
                rec.cigarstring = str(hit.q_st) + 'H' + hit.cigar_str + str(len(q.query_sequence)-hit.q_en) + 'H'
            if hit.strand == -1: 
                rec.is_reverse = True
                rec.query_sequence = mp.revcomp(q.query_sequence[hit.q_st:hit.q_en])
                rec.query_qualities = q.query_qualities[hit.q_st:hit.q_en][::-1]
                rec.cigarstring =  str(len(q.query_sequence)-hit.q_en) +'H' + hit.cigar_str + str(hit.q_st) + 'H'

            ref_olap = recs_list[0].get_overlap(hit.r_st, hit.r_en) 
            hit_refRange = hit.r_en - hit.r_st
            rlen_olap_ratio = ref_olap/hit_refRange
            if rlen_olap_ratio < 0.5:
                rec.is_qcfail = True
                recs_list[0].is_qcfail = True
            elif rec.is_supplementary and hits[0].strand == hit.strand:
                rec.is_qcfail = True
                recs_list[0].is_qcfail = True

            recs_list.append(rec)
        ## creating tags
        cigar_props = list(set(map(lambda x:x[1], hit.cigar)))
        cig_dict = dict(zip(cigar_props, [0]*len(cigar_props)))
        de_denom = 0
        for i in hit.cigar:
            cig_dict[i[1]] += i[0]
            if i[1] == 0:
                de_denom += i[0]
            else:
                de_denom += 1
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

        sa_tags.append(f"{hit.ctg},{hit.r_st},{strand_char},{str(hit.q_st) + 'S' + str(cig_dict.get(0)) + 'M' + str(cig_dict.get(2)) + 'D' + str(len(q.query_sequence) - hit.q_en) +'S'},{hit.mapq},{hit.NM}")
        tags.append([("NM", hit.NM), ("nn", (hit.blen - hit.mlen - hit.NM), "i"),("tp", tp, "A"),("de", rounded_de, "f"),("MD", hit.MD), ("cs", hit.cs)])
    
    return recs_list, tags, sa_tags


def mappy_hits_fastq_fmt(q: dict, hits: list, header_dict: dict):
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
        cigar_props = list(set(map(lambda x:x[1], hit.cigar)))
        cig_dict = dict(zip(cigar_props, [0]*len(cigar_props)))
        de_denom = 0
        for d in hit.cigar:
            cig_dict[d[1]] += d[0]
            if d[1] == 0:
                de_denom += d[0]
            else:
                de_denom += 1
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



def bam_mappy_al(bam: str, cpus: int, tech: str, ref: str,  outf: str):

    print(f"Processing {bam}")
    aligner = mappy_al_ref(ref, tech, cpus)
    
    sq_head = pysam.AlignmentFile(bam, 'rb', check_sq=False).header.to_dict() | \
    pysam.AlignmentHeader(
    ).from_references(
        reference_names=pysam.FastaFile(ref).references, 
        reference_lengths=pysam.FastaFile(ref).lengths
        ).to_dict() 
    sq_head['PG'].append({'ID': 'aligner', 'PN': 'mappy', 'VN': mp.__version__, 'DS': 'minimap2 alignment'})

    with pysam.AlignmentFile(outf, 'wb', header=sq_head) as obam:  
        for i in pysam.AlignmentFile(bam, 'rb', check_sq=False):
            hits = list(aligner.map(i.query_sequence, cs=True, MD=True))

            if not hits:
                i.is_unmapped = True
                obam.write(i)
                continue

            recs_list, tags, sa_tags = mappy_hits_bam_fmt(i, hits, sq_head)
                
            for read_alignment in recs_list:
                read_sa_tag = sa_tags.copy()
                read_sa_tag.pop(recs_list.index(read_alignment))
                if not read_sa_tag:
                    read_alignment.tags =  i.get_tags() + tags[recs_list.index(read_alignment)]
                    obam.write(read_alignment)
                    continue

                read_alignment.tags =  i.get_tags() + tags[recs_list.index(read_alignment)] + [("SA", ';'.join(read_sa_tag))]
                obam.write(read_alignment)


def fastq_mappy_al(fastq1: str, fastq2: str, cpus: int, tech: str, ref: str,  outf: str):
    print(f"Processing {fastq1}, {fastq2}")
    aligner = mappy_al_ref(ref, tech, cpus)
    
    sq_head = pysam.AlignmentHeader(
            ).from_references(
                reference_names=pysam.FastaFile(ref).references, 
                reference_lengths=pysam.FastaFile(ref).lengths
                ).to_dict()
    sq_head['PG'] = [{'ID': 'aligner', 'PN': 'mappy', 'VN': mp.__version__, 'DS': 'minimap2 alignment'}]

    with pysam.AlignmentFile(outf, 'wb', header=sq_head) as obam:
        with pysam.FastxFile(fastq1) as fq1, pysam.FastxFile(fastq2) as fq2:
            for read1, read2 in zip(fq1, fq2):
                hits = list(aligner.map(read1.sequence, seq2=read2.sequence, cs=True, MD=True))
                read_name = list(set([read1.name.split("/")[0], read2.name.split("/")[0]]))[0]

                if not hits:
                    read1.is_read1 = True
                    read2.is_read2 = True
                    for i in [read1, read2]:
                        i.is_unmapped = True
                        i.is_paired = True
                    obam.write(i)
                    continue

                rec_list = mappy_hits_fastq_fmt({"query_name": read_name, "seq1": read1.sequence, "seq2": read2.sequence, "quality1": read1.quality, "quality2": read2.quality}, hits, sq_head)

                for rec in rec_list:
                    obam.write(rec)



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
