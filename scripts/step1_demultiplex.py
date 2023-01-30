#!/usr/bin/env python3
import gzip
import click

from PROFIT_seq.seqIO import yield_fastx
from PROFIT_seq.logger import get_logger
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
    cnt = 0
    with gzip.open(outfile, 'wb') as out:
        for seq_id, seq, sep, qual in yield_fastx(infile):
            header = dict([i.split('=') for i in seq_id.split(' ')])
            ch = int(header['ch'])
            if start <= ch < end:
                out.write(f"@{seq_id}\n{seq}\n{sep}\n{qual}\n")
                cnt += 1
    LOGGER.info(f"Finished processing {cnt} reads")


if __name__ == '__main__':
    main()
