#!/usr/bin/env python3
import sys
import gzip
import click
from pathlib import Path
sys.path.append(Path(__file__).parent.parent)

from PROFIT_seq.seqIO import yield_fastx
from PROFIT_seq.logger import get_logger
from PROFIT_seq.utils import to_bytes
LOGGER = get_logger('PROFIT-seq')


@click.command()
@click.option('--infile', '-i', type=str, required=True,
              help='input gzipped fastq file.')
@click.option('--outfile', '-o', type=str, required=True,
              help='output gzipped fastq file.')
@click.option('--start', '-st', type=int, required=True,
              help='start channel number.')
@click.option('--end', '-en', type=int, required=True,
              help='end channel number.')
def main(infile, outfile, start, end):
    LOGGER.info(f"Start demultiplexing reads from channel {start}-{end}")

    cnt = 0
    with gzip.open(outfile, 'wb') as out:
        for seq_id, seq, sep, qual in yield_fastx(infile):
            header = dict([i.split('=') for i in seq_id.split(' ')[1:]])
            ch = int(header['ch'])
            if start <= ch <= end:
                out.write(to_bytes(f"@{seq_id}\n{seq}\n{sep}\n{qual}\n"))
                cnt += 1
    LOGGER.info(f"Finished processing {cnt} reads")


if __name__ == '__main__':
    main()
