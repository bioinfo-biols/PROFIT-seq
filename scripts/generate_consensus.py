#!/usr/bin/env python3
import os
import sys
import json
import click
from pathlib import Path 
from collections import defaultdict

project_dir = Path(os.path.realpath(__file__)).parent.parent
sys.path.append(str(project_dir / 'libs'))
from utils import load_fasta, get_logger
LOGGER = get_logger()


def revcomp(seq):
    """
    Convert sequence to reverse complementary
    """
    trantab = str.maketrans("ATCG", "TAGC")
    return seq.translate(trantab)[::-1]


def align_sequence(aligner, query):
    for_aln = aligner.align(query)
    rev_aln = aligner.align(revcomp(query))
    if for_aln.score > rev_aln.score:
        return for_aln
    else:
        return rev_aln


def load_adapters(adapter_fa):
    """
    Load RCA splint and 10x adapters
    """
    adapters = load_fasta(adapter_fa)
    for adapter_id in 'RCA_splint', '3prime_adapter', '5prime_adapter':
        if adapter_id not in adapters:
            sys.exit("Can not find sequence with id '{}' in adapter fasta!".format(adapter_id))
    splint = adapters['RCA_splint']
    primer3 = adapters['3prime_adapter']
    primer5 = adapters['5prime_adapter']
    return splint, primer3, primer5


def trim_adapters(ccs_fa, splint, primer3, primer5, threads):
    """
    1. Re-order consensus sequence using RCA_splint
    2. Trim 10x 5'- and 3'- adapters
    """
    from multiprocessing import Pool
    from utils import ProgressBar, grouper

    ccs_reads = []
    with open(ccs_fa, 'r') as f:
        for line in f:
            header = line.rstrip().lstrip('>')
            read_id, segments, ccs_len = header.split('\t')
            seq = f.readline().rstrip()
            ccs_reads.append([read_id, segments, int(ccs_len), seq])

    pool = Pool(threads)
    jobs = []
    for tmp in grouper(ccs_reads, 10000):
        chunk = [i for i in tmp if i is not None]
        jobs.append(pool.apply_async(trim_chunk, (chunk, splint, primer3, primer5, )))
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


def trim_chunk(chunk, splint, primer3, primer5):
    import ssw
    summary, res = [], []
    for read_id, segments, seq_len, seq in chunk:
        # if read_id != '60303d12-23f5-4643-b088-44acb1a59b2d':
        #     continue

        tmp_summary = {
            'read_id': read_id,
            'segments': segments,
            'seq_len': seq_len,
        }

        # Detect splint
        splint_alinger = ssw.Aligner(seq*2, match=2, mismatch=4, gap_open=4, gap_extend=2)
        splint_aln = align_sequence(splint_alinger, splint)
        # print(repr(splint_aln))
        
        rst, ren = splint_aln.ref_begin, splint_aln.ref_end
        if ren <= seq_len:
            ordered = seq[ren:] + seq[:rst]
            tmp_summary['splint_pos'] = '{}:{}'.format(rst, ren)
            tmp_summary['splint_type'] = 'middle'
        elif rst >= seq_len:
            ordered = seq[ren-seq_len:] + seq[:rst-seq_len]
            tmp_summary['splint_pos'] = '{}:{}'.format(rst-seq_len, ren-seq_len)
            tmp_summary['splint_type'] = 'middle'
        else:
            ordered = seq[ren-seq_len:rst]
            tmp_summary['splint_pos'] = '{}:|:{}'.format(ren-seq_len, rst)
            tmp_summary['splint_type'] = 'junction'
    
        # Filter by qc tag
        splint_cov = (splint_aln.query_end - splint_aln.query_begin) / len(splint)
        if splint_cov < 0.75:
            tmp_summary['qc_tag'] = 'no_splint'
            summary.append(tmp_summary)
            continue

        if len(ordered) < 50:
            tmp_summary['qc_tag'] = 'fragment_smaller_than_50'
            summary.append(tmp_summary)
            continue

        tmp_summary['fragment_len'] = len(ordered)

        # Trimmed adapters
        ordered_aligner = ssw.Aligner(ordered, match=1, mismatch=1, gap_open=1, gap_extend=1,)
        aln3 = align_sequence(ordered_aligner, primer3)
        aln5 = align_sequence(ordered_aligner, primer5)
        tmp_summary['primer3'] = '{}:{}'.format(aln3.ref_begin, aln3.ref_end)
        tmp_summary['primer5'] = '{}:{}'.format(aln5.ref_begin, aln5.ref_end)

        aln3_cov = (aln3.query_end - aln3.query_begin) / len(primer3)
        aln5_cov = (aln5.query_end - aln5.query_begin) / len(primer5)
        if aln3_cov < 0.6 and aln5_cov < 0.6:
            tmp_summary['qc_tag'] = 'no_both_adapters'
            summary.append(tmp_summary)
            continue
        elif aln3_cov < 0.6:
            tmp_summary['qc_tag'] = 'no_3prime_adapter'
            summary.append(tmp_summary)
            continue
        elif aln5_cov < 0.6:
            tmp_summary['qc_tag'] = 'no_5prime_adapter'
            summary.append(tmp_summary)
            continue
        else:
            pass

        # Generate consensus reads
        adapters = sorted([[aln3.ref_begin, aln3.ref_end], [aln5.ref_begin, aln5.ref_end]])
        trimmed = ordered[adapters[0][1]:adapters[1][0]]
        if len(trimmed) < 50:
            tmp_summary['qc_tag'] = 'clean_smaller_than_50'
            summary.append(tmp_summary)
            continue

        tmp_summary['qc_tag'] = 'pass'  
        tmp_summary['clean_len'] = len(trimmed)
        summary.append(tmp_summary)

        res.append([read_id, trimmed])

    return res, summary


def write_consensus(out_file, cleaned_reads, summary_file, reads_summary):
    from collections import Counter
    with open(out_file, 'w') as out:
        for read_id, trimmed in cleaned_reads:
            out.write('>{}\n{}\n'.format(read_id, trimmed))
    
    failed_reasons = []
    summary_cols = [
        'read_id', 'segments', 'seq_len', 'splint_pos', 'splint_type',
        'fragment_len', 'primer3', 'primer5', 'clean_len', 'qc_tag',
    ]
    with open(summary_file, 'w') as out:
        out.write('\t'.join(summary_cols) + '\n')
        for tmp_summary in reads_summary:
            tmp_line = [tmp_summary[i] if i in tmp_summary else 'NA' for i in summary_cols]
            out.write('\t'.join([str(x) for x in tmp_line]) + '\n')
            failed_reasons.append(tmp_summary['qc_tag'])
    return dict(Counter(failed_reasons))


@click.command()
@click.option('--ccs', '-i', type=str, required=True,
              help='CCS output from circtools.')
@click.option('--out', '-o', type=str, required=True,
              help='consensus output file.')
@click.option('--summary', '-s', type=str, required=True,
              help='qc summary file.')
@click.option('--adapter', '-r', type=str, default=project_dir / 'RCA_adapters.fa', 
              help='Adapter sequences file.')
@click.option('--threads', '-t', type=int, default=os.cpu_count(),
              help='number of threads.')            
def main(adapter, ccs, out, summary, threads):
    LOGGER.info('Processing consensus reads with {} threads'.format(threads))

    # Load adapter
    splint, primer3, primer5 = load_adapters(adapter)
    
    # Trim adapters
    cleaned_reads, reads_summary = trim_adapters(ccs, splint, primer3, primer5, threads)

    # Output
    status = write_consensus(out, cleaned_reads, summary, reads_summary)
    for i, j in status.items():
        LOGGER.info('Tag: {}, numbers: {}'.format(i, j))


if __name__ == '__main__':
    main()
