#!/usr/bin/env python3
import os
import sys
import click
import numpy as np
from pathlib import Path 
from multiprocessing import Pool, Lock
from collections import namedtuple
Adapter = namedtuple('Adapter', ['seq', 'aligner'])
LOCK = Lock()

import ssw
from spoa import poa
from pyccs import find_consensus

from logger import get_logger, ProgressBar
from seqIO import load_fastx, yield_fastx, revcomp
from utils import grouper, align_sequence
LOGGER = get_logger()


def load_adapters(adapter_fa):
    """
    Load RCA splint and 10x adapters
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
    LOGGER.error(f"{Path(__file__).name} - Error in worker process: {e.__cause__}")


def _ret_callback(ret):
    """
    Callback for updating output files in multiprocessing
    """
    # global Processed, Output
    # read_num, ret_lines = ret
    with LOCK:
        pass
        # for line in ret_lines:
        #     Output.write('\t'.join([str(x) for x in line]) + '\n')
        # Processed += read_num
        # Prog.update(100 * Processed / Total, is_force=True)


def split_and_clean(input_fa, splint, primer3, primer5, threads):
    """
    1. Re-order consensus sequence using RCA_splint
    2. Trim 10x 5'- and 3'- adapters
    """
    from multiprocessing import Pool

    pool = Pool(threads)
    jobs = []

    chunk_size = 10_000
    chunk = []
    for record in yield_fastx(input_fa):
        chunk.append(record)
        if len(chunk) >= chunk_size:
            jobs.append(pool.apply_async(worker, (chunk, splint, primer3, primer5,),
                                         callback=_ret_callback, error_callback=_err_callback, ))
            chunk = []
            break
    pool.close()

    res, summary = [], []
    n = 0
    prog = ProgressBar()
    prog.update(0)
    for job in jobs:
        tmp_res, tmp_summary = job.get()
        res += tmp_res
        summary += tmp_summary
        n += 1
        prog.update(100 * n / len(jobs))
    pool.join()
    prog.update(100)
    
    return res, summary 


def worker(chunk, splint, primer3, primer5):
    _splint = Adapter(splint, ssw.Aligner(splint, match=2, mismatch=4, gap_open=4, gap_extend=2))
    _primer3 = Adapter(primer3, ssw.Aligner(primer3, match=2, mismatch=4, gap_open=4, gap_extend=2))
    _primer5 = Adapter(primer5, ssw.Aligner(primer5, match=2, mismatch=4, gap_open=4, gap_extend=2))

    for read_id, seq, sep, qual in chunk:
        segment_str, ccs_seq = split_by_splint(seq, _splint, _primer3, _primer5)
        if ccs_seq is None:
            pass

        segment_str, ccs_seq = find_consensus(seq)
        if ccs_seq is None:
            pass



def split_by_splint(seq, splint, primer3, primer5, p_indel=0.1):
    for_aln = splint.aligner.align(seq)
    rev_aln = splint.aligner.align(revcomp(seq))
    strand_seq = revcomp(seq) if for_aln.score < rev_aln.score else seq
    splint_hits = recursive_aln_splint(strand_seq, splint, 0, True, True)
    if len(splint_hits) == 0:
        return None, None

    segments = []
    segment_str = []
    segments.append(strand_seq[:splint_hits[0][0]])
    segment_str.append(f"0:{splint_hits[0][0]}")
    if len(splint_hits) > 1:
        for l_hit, r_hit in zip(splint_hits[:-1], splint_hits[1:]):
            segments.append(strand_seq[l_hit[1]:r_hit[0]])
            segment_str.append(f"{l_hit[1]}:{r_hit[0]}")
    segments.append(strand_seq[splint_hits[-1][1]:])
    segment_str.append(f"{splint_hits[-1][1]}")

    if len(segments) > 3:
        d_estimate = np.array([len(i) for i in segments][1:-1])
        d_mean = np.mean(d_estimate)
        d_delta = np.sqrt(2.3 * (p_indel * (d_mean)))
        if np.abs(d_estimate-d_mean).max() > d_delta:
            return segment_str, None

    ccs, msa = poa(sorted(segments, key=lambda x: len(x), reverse=True))
    return segment_str, ccs


def recursive_aln_splint(strand_seq, splint, shift=0, st_align=True, en_align=True):
    if strand_seq is None:
        return []

    hit = splint.aligner.align(strand_seq)

    # Middle alignment
    #          ^-------------------^
    # =========^====================^============
    if hit.ref_end-hit.ref_begin > 0.75*len(splint.seq):
        if hit.query_begin >= 50:
            st_clip = strand_seq[:hit.query_begin]
        else:
            st_clip = None
        if hit.query_end <= len(strand_seq)-50:
            en_clip = strand_seq[hit.query_end:]
        else:
            en_clip = None
        return recursive_aln_splint(st_clip, splint, shift, st_align, False) + \
               [(hit.query_end+shift, hit.query_end+shift, hit.ref_begin, hit.ref_end)] + \
               recursive_aln_splint(en_clip, splint, shift + hit.query_end, False, en_align)

    # Head-to-head alignment
    # ----------^---------^
    #           ^=========^====================
    if st_align and hit.ref_begin >= 5 and hit.query_begin <= 5 and hit.ref_end-hit.ref_begin >= 10:
        en_clip = strand_seq[hit.query_end:]
        return [(hit.query_end+shift, hit.query_end+shift, hit.ref_begin, hit.ref_end)] + \
               recursive_aln_splint(en_clip, splint, shift+hit.query_end, False, en_align)

    # Tail-to-Tail alignment
    #                     ^----------^---------
    # ====================^==========^
    if en_align and hit.ref_begin <= 5 and hit.query_begin >= len(strand_seq)-5 and hit.ref_end-hit.ref_begin >= 10:
        st_clip = strand_seq[:hit.query_begin]
        return recursive_aln_splint(st_clip, splint, shift, st_align, False) + \
               [(hit.query_end+shift, hit.query_end+shift, hit.ref_begin, hit.ref_end)]

    # No hit
    return []


def rearrange_by_splint(strand_seq, splint):
    tmp_summary = {}
    seq_len = len(strand_seq)
    pseudo_aligner = ssw.Aligner(strand_seq*2, match=2, mismatch=4, gap_open=4, gap_extend=2)
    splint_aln = align_sequence(pseudo_aligner, splint, strandness=False)

    rst, ren = splint_aln.ref_begin, splint_aln.ref_end
    if ren <= seq_len:
        ordered = strand_seq[ren:] + strand_seq[:rst]
        tmp_summary['splint_pos'] = '{}:{}'.format(rst, ren)
        tmp_summary['splint_type'] = 'middle'
    elif rst >= seq_len:
        ordered = strand_seq[ren-seq_len:] + strand_seq[:rst-seq_len]
        tmp_summary['splint_pos'] = '{}:{}'.format(rst-seq_len, ren-seq_len)
        tmp_summary['splint_type'] = 'middle'
    else:
        ordered = strand_seq[ren-seq_len:rst]
        tmp_summary['splint_pos'] = '{}:|:{}'.format(ren-seq_len, rst)
        tmp_summary['splint_type'] = 'junction'

    # Filter by qc tag
    splint_cov = (splint_aln.query_end - splint_aln.query_begin) / len(splint)
    if splint_cov < 0.75:
        tmp_summary['qc_tag'] = 'no_splint'
        return None, tmp_summary

    if len(ordered) < 50:
        tmp_summary['qc_tag'] = 'fragment_smaller_than_50'
        return None, tmp_summary

    tmp_summary['fragment_len'] = len(ordered)
    return ordered, tmp_summary


# def trim_ccs(chunk, splint, primer3, primer5):
#     import ssw
#     summary, res = [], []
#     for read_id, segments, seq_len, seq in chunk:
#         tmp_summary = {
#             'read_id': read_id,
#             'segments': segments,
#             'seq_len': seq_len,
#         }
#
#         # Detect splint
#         splint_alinger = ssw.Aligner(seq*2, match=2, mismatch=4, gap_open=4, gap_extend=2)
#         splint_aln = align_sequence(splint_alinger, splint)
#         # print(repr(splint_aln))
#
#         rst, ren = splint_aln.ref_begin, splint_aln.ref_end
#         if ren <= seq_len:
#             ordered = seq[ren:] + seq[:rst]
#             tmp_summary['splint_pos'] = '{}:{}'.format(rst, ren)
#             tmp_summary['splint_type'] = 'middle'
#         elif rst >= seq_len:
#             ordered = seq[ren-seq_len:] + seq[:rst-seq_len]
#             tmp_summary['splint_pos'] = '{}:{}'.format(rst-seq_len, ren-seq_len)
#             tmp_summary['splint_type'] = 'middle'
#         else:
#             ordered = seq[ren-seq_len:rst]
#             tmp_summary['splint_pos'] = '{}:|:{}'.format(ren-seq_len, rst)
#             tmp_summary['splint_type'] = 'junction'
#
#         # Filter by qc tag
#         splint_cov = (splint_aln.query_end - splint_aln.query_begin) / len(splint)
#         if splint_cov < 0.75:
#             tmp_summary['qc_tag'] = 'no_splint'
#             summary.append(tmp_summary)
#             continue
#
#         if len(ordered) < 50:
#             tmp_summary['qc_tag'] = 'fragment_smaller_than_50'
#             summary.append(tmp_summary)
#             continue
#
#         tmp_summary['fragment_len'] = len(ordered)
#
#         # Trimmed adapters
#         ordered_aligner = ssw.Aligner(ordered, match=1, mismatch=1, gap_open=1, gap_extend=1,)
#         aln3 = align_sequence(ordered_aligner, primer3)
#         aln5 = align_sequence(ordered_aligner, primer5)
#         tmp_summary['primer3'] = '{}:{}'.format(aln3.ref_begin, aln3.ref_end)
#         tmp_summary['primer5'] = '{}:{}'.format(aln5.ref_begin, aln5.ref_end)
#
#         aln3_cov = (aln3.query_end - aln3.query_begin) / len(primer3)
#         aln5_cov = (aln5.query_end - aln5.query_begin) / len(primer5)
#         if aln3_cov < 0.6 and aln5_cov < 0.6:
#             tmp_summary['qc_tag'] = 'no_both_adapters'
#             summary.append(tmp_summary)
#             continue
#         elif aln3_cov < 0.6:
#             tmp_summary['qc_tag'] = 'no_3prime_adapter'
#             summary.append(tmp_summary)
#             continue
#         elif aln5_cov < 0.6:
#             tmp_summary['qc_tag'] = 'no_5prime_adapter'
#             summary.append(tmp_summary)
#             continue
#         else:
#             pass
#
#         # Generate consensus reads
#         adapters = sorted([[aln3.ref_begin, aln3.ref_end], [aln5.ref_begin, aln5.ref_end]])
#         trimmed = ordered[adapters[0][1]:adapters[1][0]]
#         if len(trimmed) < 50:
#             tmp_summary['qc_tag'] = 'clean_smaller_than_50'
#             summary.append(tmp_summary)
#             continue
#
#         tmp_summary['qc_tag'] = 'pass'
#         tmp_summary['clean_len'] = len(trimmed)
#         summary.append(tmp_summary)
#
#         res.append([read_id, trimmed])
#
#     return res, summary


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


@click.command()
@click.option('--input', '-i', type=str, required=True,
              help='CCS output from circtools.')
@click.option('--ccs', '-c', type=str, required=True,
              help='cyclic consensus reads file.')
@click.option('--non', '-n', type=str, required=True,
              help='no-cyclic consensus reads file.')
@click.option('--summary', '-s', type=str, required=True,
              help='qc summary file.')
@click.option('--adapter', '-r', type=str, default=Path(os.path.realpath(__file__)).parent.parent / 'RCA_adapters.fa',
              help='Adapter sequences file.')
@click.option('--threads', '-t', type=int, default=os.cpu_count(),
              help='number of threads.')            
def main(adapter, input, ccs, non, summary, threads):
    LOGGER.info('Processing consensus reads with {} threads'.format(threads))

    # Load adapter
    splint, primer3, primer5 = load_adapters(adapter)

    # Trim adapters
    split_and_clean(input, splint, primer3, primer5, threads)

    # Trim adapters
    # ccs_reads, non_reads, reads_summary = trim_adapters(input, splint, primer3, primer5, threads)

    # Output
    # status = write_consensus(out, cleaned_reads, summary, reads_summary)
    # for i, j in status.items():
    #     LOGGER.info('Tag: {}, numbers: {}'.format(i, j))


if __name__ == '__main__':
    main()
