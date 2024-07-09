import sys
import click
from pathlib import Path

from PROFIT_seq.logger import get_logger
Logger = get_logger("PROFIT_seq")


@click.command()
@click.option('--minknow_host', type=str, default='127.0.0.1',
              help='ip address for MinKNOW host')
@click.option('--minknow_port', type=str, default='8000',
              help='port for MinKNOW service.')
@click.option('--guppy_address', type=str, default='ipc:///tmp/.guppy/5555',
              help='address for guppy server.')
@click.option('--guppy_config', type=str, default='dna_r9.4.1_450bps_fast',
              help='guppy basecalling config.')
@click.option('--dashboard_port', type=str, default='55280',
              help='guppy basecalling config.')
@click.option('--mm_idx', type=str, required=True,
              help='Minimap2 index of reference sequences')
@click.option('--gtf', type=str, default=None,
              help='reference genome annotation')
@click.option('--debug', is_flag=True,
              help='run in debugging mode')
@click.version_option('0.2b')
def main(minknow_host, minknow_port, guppy_address, guppy_config, dashboard_port, mm_idx, gtf, debug):
    from app.server import app
    from PROFIT_seq.env import initializer

    if Path(sys.argv[0]).name != 'PROFIT_seq':
        level = 10 if debug else 20
        Logger.setLevel(level)
        for handler in Logger.handlers:
            handler.setLevel(level)

    # mm_idx = '/home/zhangjy/data/database/GENCODE/hg38.mmi'
    initializer(minknow_host, minknow_port, guppy_address, guppy_config, mm_idx, gtf)

    Logger.info("Starting PROFIT-seq dashboard")
    app.run(debug=False, host="0.0.0.0", port=dashboard_port)


if __name__ == '__main__':
    main()
