#!/usr/bin/env python3
import sys
import click
from pathlib import Path
from PROFIT_seq.logger import get_logger
from PROFIT_seq.seqIO import yield_fastx
LOGGER = get_logger("PROFIT-seq", debugging=False)


def find_polyA(seq, tail_len=30, k=10, max_A=8):
    from collections import Counter

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
@click.option('--input', '-i', type=str, required=True,
              help='input full-length consensus reads.')
@click.option('--output', '-o', type=str, required=True,
              help='output trimmed reads.')
def main(input, output):
    # Run function
    if Path(sys.argv[0]).name != 'FS-seq':
        LOGGER.setLevel(10)
        for handler in LOGGER.handlers:
            handler.setLevel(10)
    LOGGER.info(f"Start running!")

    n_reads = 0
    n_pA = 0
    with open(output, 'w') as out:
        for record in yield_fastx(input):
            n_reads += 1
            pA_st = find_polyA(record[1], 15, 10, 8)
            if pA_st is not None:
                seq = record[1][:pA_st]
                n_pA += 1
            else:
                seq = record[1]
            out.write(f">{record[0]}\n{seq}\n")

    LOGGER.info(f"Finished! Trimmed {n_pA}/{n_reads} ({100*n_pA/n_reads:.2f}%) reads with poly(A) tails!")


if __name__ == '__main__':
    main()