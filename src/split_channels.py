#!/usr/bin/env python3
import os
import sys
import click
from pathlib import Path

project_dir = Path(os.path.realpath(__file__)).parent.parent
sys.path.append(str(project_dir / 'libs'))
from utils import load_fasta, get_logger
LOGGER = get_logger()

@click.command()
@click.option('--summary', '-s', type=str, required=True,
              help='input sequencing summary.')
@click.option('--control', '-c', type=str, required=True,
              help='control file')
@click.option('--enriched', '-e', type=str, required=True,
              help='enriched file')
def main(summary, control, enriched):
    LOGGER.info("Loading sequencing summary file...")
    cnt = 0
    with open(summary, 'r') as f, open(control, 'w') as c, open(enriched, 'w') as e:
        f.readline()
        for line in f:
            content = line.rstrip().split('\t')
            read_id = content[2]
            ch = int(content[4])
            if ch % 2 == 0:
                c.write("{}\n".format(read_id))
            else:
                e.write("{}\n".format(read_id))
            cnt += 1
            if cnt % 1000000 == 0:
                LOGGER.info("Loaded {},000,000 reads...".format(cnt // 1000000))
    LOGGER.info("Finished processing {} reads...".format(cnt))


if __name__ == '__main__':
    main()
