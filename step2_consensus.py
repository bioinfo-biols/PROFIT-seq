#!/usr/bin/env python3
import os
import sys
import click
import pandas as pd
import numpy as np
from pathlib import Path
from traceback import format_exc
from multiprocessing import Lock
from collections import namedtuple
Splint_hit = namedtuple('Splint_hit', ['query_begin', 'query_end', 'ref_begin', 'ref_end', 'score'])
LOCK = Lock()

import edlib
import parasail
from pyabpoa import msa_aligner
from pyccs import find_consensus

from PROFIT_seq.logger import get_logger, ProgressBar
from PROFIT_seq.seqIO import count_fastx, load_fastx, yield_fastx, revcomp
from PROFIT_seq.utils import check_dir
LOGGER = get_logger("PROFIT-seq", debugging=False)
Columns = [
    'Read_id', "Read_type", 'Seq_len', 'Splint_hit', 'Splint_strand', 'Segment_len',
    'Consensus_split', 'Consensus_splint_strand', 'Consensus_splint_pos', 'Consensus_splint_type', 'Consensus_qc_tag',
    'Strand_len', 'Primer3', 'Primer5', 'Primer_qc_tag', 'Cleaned_len', "Blacklist", "pA", "Comments"
]
FL_Output, NonFL_output = None, None
Processed, Total = 0, 0
FL_reads, NonFL_reads = 0, 0
Prog = ProgressBar()


def load_adapters(adapter_fa):
    """
    Load RCA splint and 10x adapter sequences
    """
    adapters = load_fastx(adapter_fa)
    for adapter_id in 'RCA_splint', '3prime_adapter', '5prime_adapter':
        if adapter_id not in adapters:
            sys.exit("Can not find sequence with id '{}' in adapter fasta!".format(adapter_id))
    splint = adapters['RCA_splint']
    primer3 = adapters['3prime_adapter']
    primer5 = adapters['5prime_adapter']
    return splint, primer3, primer5


def _err_callback(e):
    """
    Callback for returning exception in multiprocessing
    """
    LOGGER.error(f"{Path(__file__).name} - Error in worker process {os.getpid()}: {e.__cause__}\n{format_exc()}")


def _ret_callback(ret):
    """
    Callback for updating output files in multiprocessing
    """
    global Processed, Output, FL_reads, NonFL_reads
    with LOCK:
        for cleaned_seq, segment_seqs, read_summary in ret:
            read_id = read_summary['Read_id'].split(' ')[0] + "|" + read_summary['UMI']
            tmp_line = [read_summary[col] if col in read_summary else "None" for col in Columns]
            tmp_line[0] = read_id
            Summary_output.write("\t".join([str(x) for x in tmp_line]) + "\n")

            if cleaned_seq is not None:
                if read_summary["Read_type"] in ["FL-splint", "FL-CCS"]:
                    FL_Output.write(f">{read_id}\n{cleaned_seq}\n")
                    FL_reads += 1
                else:
                    NonFL_output.write(f">{read_id}\n{cleaned_seq}\n")
                    NonFL_reads += 1

        # LOGGER.debug(f"Processed {n_fl}/{n_nonfl}/{n_total} ({100*n_fl/n_total:.2f}%) reads")
        Processed += len(ret)
        Prog.update(100 * Processed / Total, is_force=True)


def split_and_clean(input_fa, splint, primer3, primer5, blacklist, threads, trimA, is_denovo, is_permissive):
    """
    1. Re-order consensus sequence using RCA_splint
    2. Trim 10x 5'- and 3'- adapters
    """
    from multiprocessing import Pool
    global Total

    pool = Pool(threads)
    jobs = []

    chunk_size = 1_000
    chunk = []
    Prog.update(0)
    for record in yield_fastx(input_fa):
        chunk.append(record)
        if len(chunk) >= chunk_size:
            jobs.append(pool.apply_async(worker, (chunk, splint, primer3, primer5, blacklist, trimA, is_denovo, is_permissive),
                                         callback=_ret_callback, error_callback=_err_callback, ))
            chunk = []
    pool.close()
    pool.join()
    Prog.close()

    return 0


def worker(chunk, splint, primer3, primer5, blacklist, trimA, is_denovo=False, is_permissive=False):
    """
    Main worker process
    """
    aligner = msa_aligner()
    umi_pos = find_umi(splint, "NNNNN")
    ret = []
    for read_id, seq, sep, qual in chunk:
        if len(seq) == 0:
            continue
        # if read_id.split(" ")[0] != "551a785c-4f8e-4c1e-8e12-aaa32c072a1b":
        #     continue
        tmp_summary = {
            'Read_id': read_id,
            'Seq_len': len(seq),
            'Blacklist': read_id in blacklist,
        }

        # Split by splint sequences
        segment_str, segment_strand, umi, consensus_seq, segment_seqs, comments = split_by_splint(seq, splint, aligner, umi_pos, is_permissive)
        tmp_summary['Comments'] = comments

        # Splint information
        tmp_summary['UMI'] = umi
        if segment_str is not None:
            tmp_summary['Splint_hit'] = segment_str
            tmp_summary['Splint_strand'] = segment_strand

        # Process consensus reads
        strand_seq = seq
        if consensus_seq is not None:
            # Full-length consensus reads
            strand_seq = consensus_seq
            tmp_summary['Segment_len'] = len(consensus_seq)
            tmp_summary["Read_type"] = "FL-splint"
        else:
            if is_denovo:
                # De novo find CCS sequence
                ccs_str, ccs_seq = find_consensus(seq)
                tmp_summary["Consensus_split"] = ccs_str
                if ccs_seq is not None:
                    tmp_summary["Read_type"] = "FL-CCS"
                    umi, strand_seq, reorder_summary = rearrange_by_splint(ccs_seq, splint, umi_pos)
                    tmp_summary['UMI'] = umi
                    tmp_summary.update(reorder_summary)
            else:
                # Get non-full length sequence
                tmp_summary["Read_type"] = "Non-FL"
                if umi is None:
                    umi = get_umi(splint, seq, umi_pos)
                    tmp_summary['UMI'] = umi

        # Trim adapters
        if strand_seq:
            tmp_summary["Strand_len"] = len(strand_seq)
            cleaned_seq, trim_summary = trim_strand_primer(strand_seq, primer3, primer5)
            tmp_summary.update(trim_summary)

            if cleaned_seq is not None:
                pA_st = find_polyA(cleaned_seq, 15, 10, 8)
                if pA_st is None:
                    tmp_summary['pA'] = False
                else:
                    tmp_summary['pA'] = True
                    if trimA:
                        cleaned_seq = cleaned_seq[:pA_st]
        else:
            tmp_summary['Strand_len'] = 0
            cleaned_seq = None

        if cleaned_seq is not None and len(cleaned_seq) == 0:
            cleaned_seq = None

        ret.append((cleaned_seq, segment_seqs, tmp_summary))

    return ret


def ssw_parasail(ref, query, match=2, mismatch=-4, gap_open=4, gap_extend=2):
    """
    SSW alignment of (read, adapter)
    """
    score_size = 1
    matrix = parasail.matrix_create("ACGT", match, mismatch)
    profile = parasail.ssw_init(ref, matrix, score_size)
    hit = parasail.ssw_profile(profile, query, gap_open, gap_extend)
    if hit is None:
        return Splint_hit(None, None, None, None, 0)

    return Splint_hit(hit.read_begin1, hit.read_end1, hit.ref_begin1, hit.ref_end1, hit.score1)


def aln_sequence(ref, query, strandness=True, **kwargs):
    """Align sequence and get the correct strand"""
    for_aln = ssw_parasail(ref, query, **kwargs)
    if strandness:
        return for_aln, "+"

    rev_aln = ssw_parasail(ref, revcomp(query), **kwargs)
    if for_aln.score > rev_aln.score:
        return for_aln, "+"
    else:
        return rev_aln, "-"


def split_by_splint(seq, splint, aligner, umi_pos, p_indel=0.1, permissive=False):
    """
    FInd splint sequence
    """
    if len(seq) == 0:
        return None, None, None, None, None, []
    for_aln = ssw_parasail(seq, splint)
    rev_aln = ssw_parasail(seq, revcomp(splint))

    strand, strand_seq = "+", seq
    if for_aln.score < rev_aln.score:
        strand = "-"
        strand_seq = revcomp(seq)

    comments = []

    # Get splint hit
    splint_hits = recursive_aln_splint(strand_seq, splint, 0, True, True)
    if len(splint_hits) == 0:
        return None, None, None, None, None, comments

    # Get Splint UMI sequence
    splint_seqs = [strand_seq[i:j] for i, j, _, _ in splint_hits]
    fl_splints = [i for i in splint_seqs if len(i) >= .75 * len(splint)]
    if len(fl_splints) > 1:
        splint_msa = aligner.msa(fl_splints, out_cons=True, out_msa=False)
        splint_ccs = splint_msa.cons_seq[0]
    elif len(fl_splints) == 1:
        splint_ccs = fl_splints[0]
    else:
        splint_ccs = sorted(splint_seqs, key=lambda x: len(x), reverse=True)[0]
    splint_umi = get_umi(splint, splint_ccs, umi_pos)

    # Get all sub reads
    segments, segment_str = [], []

    # Head
    is_head_splint = True
    if splint_hits[0][0] > 50:
        is_head_splint = False
        segments.append(strand_seq[:splint_hits[0][0]])
        segment_str.append(f"0:{splint_hits[0][0]}")

    # Middle
    if len(splint_hits) > 1:
        for l_hit, r_hit in zip(splint_hits[:-1], splint_hits[1:]):
            if r_hit[0]-l_hit[1] > 50:
                segments.append(strand_seq[l_hit[1]:r_hit[0]])
                segment_str.append(f"{l_hit[1]}:{r_hit[0]}")

    # Tail
    is_tail_splint = True
    if len(seq) - splint_hits[-1][1] > 50:
        is_tail_splint = False
        segments.append(strand_seq[splint_hits[-1][1]:])
        segment_str.append(f"{splint_hits[-1][1]}:{len(seq)}")

    # Full-length cDNA segments
    fl_st = 0 if is_head_splint else 1
    fl_en = len(segments) if is_tail_splint else len(segments)-1
    fl_segments = segments[fl_st:fl_en]

    # Filter based on edit distance
    ccs, non_fl = None, None
    if len(fl_segments) == 1:
        # splint - cDNA - splint
        ccs = fl_segments[0]
    elif len(fl_segments) > 1:
        # Multiple full-length cDNA
        msa = aligner.msa(fl_segments, out_cons=True, out_msa=False)
        ccs = msa.cons_seq[0]

        # Filter based on segment length / distance threshold
        d_estimate = np.array([len(i) for i in fl_segments])
        d_mean = np.mean(d_estimate)
        # d_delta = 2.3 * np.sqrt(p_indel * d_mean)
        d_delta = 0.2 * d_mean
        d = [edlib.align(i, ccs)['editDistance'] for i in fl_segments]
        comments = [np.abs(d_estimate-d_mean).max()/d_mean, np.max(d) / len(ccs)]
        if not permissive and ((np.abs(d_estimate-d_mean).max() > d_delta) or (np.max(d)/len(ccs) > 0.5)):
            ccs = None
    elif len(segments) == 2:
        # cDNA - splint - cDNA
        # Assemble full length sequence based on this rule:
        # head: ----------^---------^
        #                 ^=========^==================== tail
        hit = ssw_parasail(segments[0], segments[1])
        if hit.score != 0:
            head_st, head_en = hit.query_begin, hit.query_end
            tail_st, tail_en = hit.ref_begin, hit.ref_end
            if head_en >= len(segments[0]) - 10 and tail_st <= 10 and min(head_en-head_st, tail_en-tail_st) >= 50:
                ccs = segments[0][:head_st] + segments[1][tail_st:]

    return segment_str, strand, splint_umi, ccs, segments, comments


def recursive_aln_splint(strand_seq, splint, shift=0, st_align=True, en_align=True):
    """
    Align splint sequences
    """
    if strand_seq is None or len(strand_seq) == 0:
        return []

    hit = ssw_parasail(strand_seq, splint)
    if hit.score != 0:
        # Middle alignment
        #          ^-------------------^
        # =========^====================^============
        if hit.ref_end-hit.ref_begin >= 0.75*len(splint):
            if hit.query_begin >= 50:
                st_clip = strand_seq[:hit.query_begin]
            else:
                st_clip = None
            if hit.query_end <= len(strand_seq)-50:
                en_clip = strand_seq[hit.query_end:]
            else:
                en_clip = None
            return recursive_aln_splint(st_clip, splint, shift, st_align, False) + \
                   [(hit.query_begin+shift, hit.query_end+shift, hit.ref_begin, hit.ref_end)] + \
                   recursive_aln_splint(en_clip, splint, shift+hit.query_end, False, en_align)

        # Head-to-head alignment
        # ----------^---------^
        #           ^=========^====================
        if st_align and hit.ref_end >= len(splint) - 5 and hit.query_begin <= 5 and hit.ref_end-hit.ref_begin >= 20:
            en_clip = strand_seq[hit.query_end:]
            return [(hit.query_begin+shift, hit.query_end+shift, hit.ref_begin, hit.ref_end)] + \
                   recursive_aln_splint(en_clip, splint, shift+hit.query_end, False, en_align)

        # Tail-to-Tail alignment
        #                     ^----------^---------
        # ====================^==========^
        if en_align and hit.ref_begin <= 5 and hit.query_end >= len(strand_seq)-5 and hit.ref_end-hit.ref_begin >= 20:
            st_clip = strand_seq[:hit.query_begin]
            return recursive_aln_splint(st_clip, splint, shift, st_align, False) + \
                   [(hit.query_begin+shift, hit.query_end+shift, hit.ref_begin, hit.ref_end)]

    # No hit
    return []


def find_umi(template, umi="NNNNN", shift=0):
    """find umi Position in splint template"""
    shift = template.find(umi, shift+len(umi))
    if shift == -1:
        return []
    return [(shift, shift+len(umi))] + find_umi(template, umi, shift)


def decode(b):
    """
    Decode parasail cigar output
    """
    __BAM_CIGAR_STR = 'MIDNSHP=X8'

    def _decode(x):
        l = int(str(x >> 4))
        try:
            c = __BAM_CIGAR_STR[x & 0xf]
        except:
            c = 'M'
        return (l, c)

    return [_decode(i) for i in b]


def get_umi(splint_template, splint_ccs, umi_pos, match=2, mismatch=-4, gap_open=4, gap_extend=2):
    """
    Find UMI sequence in consensus splint
    """
    score_size = 1
    matrix = parasail.matrix_create("ACGT", match, mismatch)
    profile = parasail.ssw_init(splint_template, matrix, score_size)
    hit = parasail.ssw_profile(profile, splint_ccs, gap_open, gap_extend)

    template_st, template_en = hit.read_begin1, hit.read_end1
    ccs_st, ccs_en = hit.ref_begin1, hit.ref_end1

    template_to_ccs = [-1 for _ in range(template_st)]
    for l, op in decode(hit.cigar):
        if op == "=" or op == "M" or op == "X":  # Consume both
            template_to_ccs += list(range(ccs_st, ccs_st + l))
            ccs_st += l
        elif op == "D" or op == "N":  # Consume query
            ccs_st += l
        elif op == "I" or op == "S":  # Consume reference
            template_to_ccs += [-1 for _ in range(l)]
        elif op == "H" or op == "P":  # Consume None
            continue

    template_to_ccs = np.array(template_to_ccs, dtype=np.int16)
    umi_seq = []
    for st, en in umi_pos:
        _umi = ""
        for x in range(st, en):
            if x >= len(template_to_ccs):
                continue
            if template_to_ccs[x] == -1:
                continue
            _umi += splint_ccs[template_to_ccs[x]]
        umi_seq.append(_umi)
    umi_seq = ';'.join(umi_seq)
    return umi_seq


def rearrange_by_splint(strand_seq, splint, umi_pos):
    """
    Re arrange consensus sequence by splint position
    """
    tmp_summary = {}
    seq_len = len(strand_seq)
    splint_aln, splint_strand = aln_sequence(strand_seq*2, splint, strandness=False)

    rst, ren = splint_aln.query_begin, splint_aln.query_end
    umi = get_umi(splint, (strand_seq*2)[rst:ren], umi_pos)
    if ren <= seq_len:
        ordered = strand_seq[ren:] + strand_seq[:rst]
        tmp_summary['Consensus_splint_pos'] = '{}:{}'.format(rst, ren)
        tmp_summary['Consensus_splint_type'] = 'middle'
    elif rst >= seq_len:
        ordered = strand_seq[ren-seq_len:] + strand_seq[:rst-seq_len]
        tmp_summary['Consensus_splint_pos'] = '{}:{}'.format(rst-seq_len, ren-seq_len)
        tmp_summary['Consensus_splint_type'] = 'middle'
    else:
        ordered = strand_seq[ren-seq_len:rst]
        tmp_summary['Consensus_splint_pos'] = '{}:|:{}'.format(ren-seq_len, rst)
        tmp_summary['Consensus_splint_type'] = 'junction'
    tmp_summary['Consensus_splint_strand'] = splint_strand

    # Filter by qc tag
    splint_cov = (splint_aln.ref_end - splint_aln.ref_begin) / len(splint)
    if splint_cov < 0.75:
        tmp_summary['Consensus_qc_tag'] = 'no_splint'
        return umi, None, tmp_summary

    if len(ordered) < 50:
        tmp_summary['Consensus_qc_tag'] = 'fragment_smaller_than_50'
        return umi, None, tmp_summary

    return umi, ordered, tmp_summary


def trim_strand_primer(strand_seq, primer3, primer5):
    """
    Detect 3' and 5' primer
    """
    tmp_summary = {}
    # Trimmed adapters
    # strand_aligner = ssw.Aligner(strand_seq, match=1, mismatch=1, gap_open=1, gap_extend=1,)
    aln3, _ = aln_sequence(strand_seq, revcomp(primer3), strandness=True)
    aln5, _ = aln_sequence(strand_seq, primer5, strandness=True)

    tmp_summary['Primer3'] = '{}:{}'.format(aln3.query_begin, aln3.query_end)
    tmp_summary['Primer5'] = '{}:{}'.format(aln5.query_begin, aln5.query_end)

    aln3_cov = (aln3.ref_end - aln3.ref_begin) / len(primer3)
    aln5_cov = (aln5.ref_end - aln5.ref_begin) / len(primer5)
    if aln3_cov < 0.6 and aln5_cov < 0.6:
        tmp_summary['Primer_qc_tag'] = 'no_both_adapters'
        return None, tmp_summary
    elif aln3_cov < 0.6:
        tmp_summary['Primer_qc_tag'] = 'no_3prime_adapter'
        return None, tmp_summary
    elif aln5_cov < 0.6:
        tmp_summary['Primer_qc_tag'] = 'no_5prime_adapter'
        return None, tmp_summary
    else:
        pass

    # Generate consensus reads
    adapters = sorted([[aln3.query_begin, aln3.query_end], [aln5.query_begin, aln5.query_end]])
    trimmed = strand_seq[adapters[0][1]:adapters[1][0]]
    if len(trimmed) < 50:
        tmp_summary['Primer_qc_tag'] = 'clean_smaller_than_50'
        return None, tmp_summary

    tmp_summary['Primer_qc_tag'] = 'pass'
    tmp_summary['Cleaned_len'] = len(trimmed)

    return trimmed, tmp_summary


# def write_consensus(out_file, cleaned_reads, summary_file, reads_summary):
#     from collections import Counter
#     with open(out_file, 'w') as out:
#         for read_id, trimmed in cleaned_reads:
#             out.write('>{}\n{}\n'.format(read_id, trimmed))
#
#     failed_reasons = []
#     summary_cols = [
#         'read_id', 'segments', 'seq_len', 'splint_pos', 'splint_type',
#         'fragment_len', 'primer3', 'primer5', 'clean_len', 'qc_tag',
#     ]
#     with open(summary_file, 'w') as out:
#         out.write('\t'.join(summary_cols) + '\n')
#         for tmp_summary in reads_summary:
#             tmp_line = [tmp_summary[i] if i in tmp_summary else 'NA' for i in summary_cols]
#             out.write('\t'.join([str(x) for x in tmp_line]) + '\n')
#             failed_reasons.append(tmp_summary['qc_tag'])
#     return dict(Counter(failed_reasons))


def load_sequencing_summary(summary_file):
    """
    Load blacklist of unblock reads from sequencing summary
    """
    if summary_file is None:
        return {}
    df = pd.read_csv(summary_file, sep="\t", index_col=2, low_memory=False)
    blacklist = df.index[df['end_reason'].map(lambda x: "unblock_mux_change" in x)]
    return set(blacklist)


def find_polyA(seq, tail_len=30, k=10, max_A=8):
    """
    Find polyA position in input sequence
    """
    from collections import Counter
    if len(seq) == 0:
        return None

    pos_A = []
    scan_en = max(len(seq) - k, 0)
    scan_st = max(len(seq) - tail_len, 0)
    for i in range(scan_st, scan_en):
        if Counter(seq[i:i + k])["A"] >= max_A:
            pos_A.append(i)

    x = scan_st
    while Counter(seq[x:x + k])["A"] >= max_A:
        pos_A.append(x)
        x -= 1

    if len(pos_A) == 0:
        return None

    st_A = min(pos_A)
    while seq[st_A] != "A":
        st_A += 1

    return st_A


@click.command()
@click.option('--input', '-i', type=Path, required=True,
              help='input porechop trimmed fastq.')
@click.option('--summary', '-s', type=Path, required=False,
              help='input sequencing summary generate by MinKNOW.')
@click.option('--adapter', '-r', type=Path, default=Path(os.path.realpath(__file__)).parent / 'RCA_adapters.fa',
              help='Adapter sequences file. Defaults to embedded splint adapter sequences.')
@click.option('--outdir', '-o', type=Path, required=True,
              help='output directory name.')
@click.option('--prefix', '-p', type=str, required=True,
              help='output prefix name.')
@click.option('--threads', '-t', type=int, default=os.cpu_count(),
              help='number of threads. Defaults to number of cpu cores.')
@click.option('--trimA', is_flag=True, help="trim 3' poly(A) sequences")
@click.option('--denovo', is_flag=True, help="use denovo consensus mode to get better sensitivity")
@click.option('--permissive', is_flag=True, help="run in permissive mode, generates more consensus reads but do not filter potential RCA artifacts")
def main(input, summary, adapter, outdir, prefix, threads, trima, denovo, permissive):
    global Total, FL_Output, NonFL_output, Summary_output

    # Run function
    if Path(sys.argv[0]).name != 'PROFIT_seq':
        LOGGER.setLevel(10)
        for handler in LOGGER.handlers:
            handler.setLevel(10)
    LOGGER.info('Processing consensus reads with {} threads'.format(threads))

    LOGGER.info('Loading sequencing summary')
    blacklist = load_sequencing_summary(summary)

    # Load adapter
    splint, primer3, primer5 = load_adapters(adapter)

    # Open output files
    outdir = check_dir(outdir)
    fl_fa = outdir / f"{prefix}.fl.fa"
    partial_fa = outdir / f"{prefix}.recovered.fa"
    summary_txt = outdir / f"{prefix}.consensus_summary.txt"

    FL_Output = open(fl_fa, 'w')
    NonFL_output = open(partial_fa, 'w')
    Summary_output = open(summary_txt, 'w')
    Summary_output.write('\t'.join(Columns) + '\n')

    # Count total reads number
    Total = count_fastx(input)
    LOGGER.info(f"Total {Total} reads")

    # Trim adapters
    status = split_and_clean(input, splint, primer3, primer5, blacklist, threads, trima, denovo, permissive)

    LOGGER.info(f"Processed {FL_reads}/{NonFL_reads}/{Total} ({100*(FL_reads+NonFL_reads)/Total:.2f}% cleaned) reads")
    LOGGER.info(f"Finished!")
    
    # Cleaning up
    FL_Output.close()
    NonFL_output.close()
    Summary_output.close()


if __name__ == '__main__':
    main()
