#!/usr/bin/env python3
import os
import sys
import click
import numpy as np
from pathlib import Path 
from multiprocessing import Lock
from collections import namedtuple
Adapter = namedtuple('Adapter', ['seq', 'aligner'])
LOCK = Lock()

import ssw
from spoa import poa
from pyccs import find_consensus

from PROFIT_seq.logger import get_logger, ProgressBar
from PROFIT_seq.seqIO import count_fastx, load_fastx, yield_fastx, revcomp
from PROFIT_seq.utils import align_sequence
LOGGER = get_logger("FS-seq", debugging=False)
Columns = [
    'Read_id', "Read_type", 'Seq_len', 'Splint_hit', 'Segment_len',
    'Consensus_split', 'Consensus_splint_pos', 'Consensus_splint_type', 'Consensus_qc_tag', 'Strand_len',
    'Primer3', 'Primer5', 'Primer_qc_tag', 'Cleaned_len'
]
FL_Output = None
NonFL_output = None
Summary_output = None
Processed, Total = 0, 0
FL_reads, NonFL_reads = 0, 0
Prog = ProgressBar()


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
    global Processed, Output, FL_reads, NonFL_reads
    with LOCK:
        for cleaned_seq, read_summary in ret:
            read_id = read_summary['Read_id'].split(' ')[0]
            tmp_line = [read_summary[col] if col in read_summary else "-" for col in Columns]
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


def split_and_clean(input_fa, splint, primer3, primer5, threads):
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
            jobs.append(pool.apply_async(worker, (chunk, splint, primer3, primer5,),
                                         callback=_ret_callback, error_callback=_err_callback, ))
            chunk = []
    pool.close()
    pool.join()
    Prog.close()

    return 0


def worker(chunk, splint, primer3, primer5):
    _splint = Adapter(splint, ssw.Aligner(splint, match=2, mismatch=4, gap_open=4, gap_extend=2))
    _primer3 = Adapter(primer3, ssw.Aligner(primer3, match=2, mismatch=4, gap_open=4, gap_extend=2))
    _primer5 = Adapter(primer5, ssw.Aligner(primer5, match=2, mismatch=4, gap_open=4, gap_extend=2))

    ret = []
    for read_id, seq, sep, qual in chunk:
        tmp_summary = {
            'Read_id': read_id,
            'Seq_len': len(seq),
        }

        # Split by splint sequences
        segment_str, consensus_seq = split_by_splint(seq, _splint)
        if segment_str is None:
            tmp_summary['Splint_hit'] = "No splint hits"
        else:
            tmp_summary['Splint_hit'] = segment_str

        if consensus_seq is None:
            tmp_summary['Segment_len'] = "Not detected"
        else:
            tmp_summary['Segment_len'] = len(consensus_seq)

        if consensus_seq is not None:
            strand_seq = consensus_seq
            tmp_summary["Read_type"] = "FL-splint"
        else:
            # Split by de novo CCS
            ccs_str, ccs_seq = find_consensus(seq)
            if ccs_str is not None:
                tmp_summary["Consensus_split"] = ccs_str
            else:
                tmp_summary["Consensus_split"] = "No CCS detected"

            # Fall back to full-length sequence
            if ccs_seq is not None:
                tmp_summary["Read_type"] = "FL-CCS"
                strand_seq, reorder_summary = rearrange_by_splint(ccs_seq, splint)
                tmp_summary.update(reorder_summary)
            else:
                tmp_summary["Read_type"] = "Non-FL"
                strand_seq = seq

        # Trim adapters
        if strand_seq:
            tmp_summary["Strand_len"] = len(strand_seq)
            cleaned_seq, trim_summary = trim_strand_primer(strand_seq, primer3, primer5)
            tmp_summary.update(trim_summary)
        else:
            tmp_summary['Strand_len'] = 0
            cleaned_seq = None

        ret.append((cleaned_seq, tmp_summary))

    return ret


def split_by_splint(seq, splint, p_indel=0.1):
    for_aln = splint.aligner.align(seq)
    rev_aln = splint.aligner.align(revcomp(seq))
    strand_seq = revcomp(seq) if for_aln.score < rev_aln.score else seq
    splint_hits = recursive_aln_splint(strand_seq, splint, 0, True, True)
    if len(splint_hits) == 0:
        return None, None

    segments = []
    segment_str = []
    if splint_hits[0][0] > 50:
        segments.append(strand_seq[:splint_hits[0][0]])
        segment_str.append(f"0:{splint_hits[0][0]}")

    if len(splint_hits) > 1:
        for l_hit, r_hit in zip(splint_hits[:-1], splint_hits[1:]):
            segments.append(strand_seq[l_hit[1]:r_hit[0]])
            segment_str.append(f"{l_hit[1]}:{r_hit[0]}")

    if len(seq) - splint_hits[-1][1] > 50:
        segments.append(strand_seq[splint_hits[-1][1]:])
        segment_str.append(f"{splint_hits[-1][1]}:{len(seq)}")

    if len(segments) > 3:
        d_estimate = np.array([len(i) for i in segments][1:-1])
        d_mean = np.mean(d_estimate)
        d_delta = np.sqrt(2.3 * (p_indel * (d_mean)))
        if np.abs(d_estimate-d_mean).max() > d_delta:
            return segment_str, None

    ccs, _ = poa(sorted(segments, key=lambda x: len(x), reverse=True))

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
               [(hit.query_begin+shift, hit.query_end+shift, hit.ref_begin, hit.ref_end)] + \
               recursive_aln_splint(en_clip, splint, shift+hit.query_end, False, en_align)

    # Head-to-head alignment
    # ----------^---------^
    #           ^=========^====================
    if st_align and hit.ref_begin >= 5 and hit.query_begin <= 5 and hit.ref_end-hit.ref_begin >= 10:
        en_clip = strand_seq[hit.query_end:]
        return [(hit.query_begin+shift, hit.query_end+shift, hit.ref_begin, hit.ref_end)] + \
               recursive_aln_splint(en_clip, splint, shift+hit.query_end, False, en_align)

    # Tail-to-Tail alignment
    #                     ^----------^---------
    # ====================^==========^
    if en_align and hit.ref_begin <= 5 and hit.query_begin >= len(strand_seq)-5 and hit.ref_end-hit.ref_begin >= 10:
        st_clip = strand_seq[:hit.query_begin]
        return recursive_aln_splint(st_clip, splint, shift, st_align, False) + \
               [(hit.query_begin+shift, hit.query_end+shift, hit.ref_begin, hit.ref_end)]

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

    # Filter by qc tag
    splint_cov = (splint_aln.query_end - splint_aln.query_begin) / len(splint)
    if splint_cov < 0.75:
        tmp_summary['Consensus_qc_tag'] = 'no_splint'
        return None, tmp_summary

    if len(ordered) < 50:
        tmp_summary['Consensus_qc_tag'] = 'fragment_smaller_than_50'
        return None, tmp_summary

    return ordered, tmp_summary


def trim_strand_primer(strand_seq, primer3, primer5):
    tmp_summary = {}
    # Trimmed adapters
    strand_aligner = ssw.Aligner(strand_seq, match=1, mismatch=1, gap_open=1, gap_extend=1,)
    aln3 = align_sequence(strand_aligner, revcomp(primer3), strandness=True)
    aln5 = align_sequence(strand_aligner, primer5, strandness=True)

    tmp_summary['Primer3'] = '{}:{}'.format(aln3.ref_begin, aln3.ref_end)
    tmp_summary['Primer5'] = '{}:{}'.format(aln5.ref_begin, aln5.ref_end)

    aln3_cov = (aln3.query_end - aln3.query_begin) / len(primer3)
    aln5_cov = (aln5.query_end - aln5.query_begin) / len(primer5)
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
    adapters = sorted([[aln3.ref_begin, aln3.ref_end], [aln5.ref_begin, aln5.ref_end]])
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


@click.command()
@click.option('--input', '-i', type=str, required=True,
              help='input trimmed fastq.')
@click.option('--ccs', '-c', type=str, required=True,
              help='output full-length consensus reads.')
@click.option('--non', '-n', type=str, required=True,
              help='output non full-length reads.')
@click.option('--summary', '-s', type=str, required=True,
              help='consensus calling summary file.')
@click.option('--adapter', '-r', type=str, default=Path(os.path.realpath(__file__)).parent / 'RCA_adapters.fa',
              help='Adapter sequences file.')
@click.option('--threads', '-t', type=int, default=os.cpu_count(),
              help='number of threads.')            
def main(adapter, input, ccs, non, summary, threads):
    global Total, FL_Output, NonFL_output, Summary_output

    # Run function
    if Path(sys.argv[0]).name != 'PROFIT_seq':
        LOGGER.setLevel(10)
        for handler in LOGGER.handlers:
            handler.setLevel(10)

    LOGGER.info('Processing consensus reads with {} threads'.format(threads))

    # Load adapter
    splint, primer3, primer5 = load_adapters(adapter)

    FL_Output = open(ccs, 'w')
    NonFL_output = open(non, 'w')
    Summary_output = open(summary, 'w')

    Total = count_fastx(input)
    LOGGER.info(f"Total {Total} reads")

    # Trim adapters
    status = split_and_clean(input, splint, primer3, primer5, threads)

    FL_Output.close()
    NonFL_output.close()
    Summary_output.close()

    LOGGER.info(f"Processed {FL_reads}/{NonFL_reads}/{Total} ({100*(FL_reads+NonFL_reads)/Total:.2f}% cleaned) reads")
    LOGGER.info(f"Finished!")


if __name__ == '__main__':
    main()
