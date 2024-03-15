"""Tools of PROFIT-seq
This module contain tools functions of PROFIT-seq, the core functions are
inspired by Readfish(https://github.com/LooseLab/readfish) and
read_until_api (https://github.com/nanoporetech/read_until_api)

Todo:
    * Add transcript coverage determine for RNA-seq

"""
import os
import sys
import time
import toml
import hashlib
import numpy as np
import pandas as pd
from pathlib import Path
from enum import IntEnum
from logging import getLogger
from itertools import zip_longest
from collections import defaultdict, Counter

from minknow_api import Connection
from pyguppy_client_lib.pyclient import PyGuppyClient
from pyguppy_client_lib.helper_functions import package_read
Logger = getLogger('PROFIT_seq')

from PROFIT_seq.parser import yield_gff


class Severity(IntEnum):
    """Severity level of MinKNOW API"""
    INFO = 1
    WARN = 2
    ERROR = 3


def send_message(rpc_connection: Connection, message: str, severity: Severity):
    """Send a message to MinKNOW
    From https://github.com/LooseLab/readfish/ru/utils.py

    Parameters
    ----------
    rpc_connection
        An instance of the rpc.Connection
    message : str
        The message to send
    severity : int
        The severity to use for the message: 1=info, 2=warning, 3=error
    Returns
    -------
    None
    """
    rpc_connection.log.send_user_message(severity=severity, user_message=message)


def basecall(guppy_client: PyGuppyClient, reads: list, dtype: "np.dtype", daq_values: dict):
    """Generator that sends and receives data from guppy
    :param guppy_client: pyguppy_client_lib.pyclient.PyGuppyClient
    :param reads: List of reads from read_until
    :type reads: Iterable
    :param dtype:
    :param daq_values:
    :returns:
        - read_info (:py:class:`tuple`) - channel (int), read number (int)
        - read_data (:py:class:`dict`) - Data returned from Guppy
    :rtype: Iterator[tuple[tuple, dict]]
    """
    hold = {}
    missing_count = 0

    with guppy_client:
        for channel, read in reads:
            hold[read.id] = (channel, read.number)
            t_0 = time.time()
            success = guppy_client.pass_read(
                package_read(
                    read_id=read.id,
                    raw_data=np.frombuffer(read.raw_data, dtype),
                    daq_offset=daq_values[channel].offset,
                    daq_scaling=daq_values[channel].scaling,
                )
            )
            if not success:
                Logger.warning("Skipped a read: %s", read.id)
                hold.pop(read.id)
                continue
            else:
                missing_count += 1

            sleep_time = guppy_client.throttle - t_0
            if sleep_time > 0:
                time.sleep(sleep_time)

        while missing_count:
            results = guppy_client.get_completed_reads()
            missing_count -= len(results)

            if not results:
                time.sleep(guppy_client.throttle)
                continue

            yield from iter(
                [(hold[read[0]["metadata"]["read_id"]], read) for read in results]
            )


def convert_time(duration: int):
    """Convert time to formatted string"""
    daytime = 86400
    hourtime = 3600
    mintime = 60

    tmp_time = duration

    n_day = duration // daytime
    tmp_time -= n_day * daytime

    n_hour = tmp_time // hourtime
    tmp_time -= n_hour * hourtime

    n_min = tmp_time // mintime
    tmp_time -= n_min * mintime

    if n_day > 0:
        return f'{n_day:.0f}d{n_hour:02.0f}h{n_min:02.0f}m{tmp_time:02.0f}s'
    elif n_hour > 0:
        return f'{n_hour:.0f}h{n_min:02.0f}m{tmp_time:02.0f}s'
    else:
        return f'{n_min:.0f}m{tmp_time:02.0f}s'


def load_toml(infile):
    """Load jobs from TOML files"""
    lines = []
    for line in infile:
        lines.append(to_str(line).strip())
    toml_str = '\n'.join(lines)
    dat = toml.loads(toml_str)
    assert "jobs" in dat

    for job in dat['jobs']:
        actions = []
        job['time'] = tuple(job['time'])
        job['ch'] = tuple(job['ch'])
        region_idx = defaultdict(list)
        for tgt in job['target']:
            if tgt['region'] not in ["multi", "mapped", "miss", "unmapped", "all"]:
                chrom, st, en = tgt['region']
                region_idx[chrom].append((int(st), int(en), tgt['action']))
            else:
                if tgt['region'] in region_idx:
                    Logger.error(f"Conflict job for {job}")
                    sys.exit()
                region_idx[tgt['region']] = tgt['action']
            actions.append(tgt['action'])

        counter = Counter(actions)
        n_max = counter.most_common(n=1)[0][1]
        a_max = [i for i, j in counter.items() if j == n_max]
        for x in ['stop_receiving', 'unblock', 'balance', 'wait']:
            if x in a_max:
                job['group'] = x
                break
        assert 'group' in job

        # Sort region
        for contig in region_idx:
            if contig not in ["multi", "mapped", "miss", "unmapped", "all"]:
                region_idx[contig] = sorted(region_idx[contig])

        job['index'] = region_idx

    return dat['jobs']


def load_prim_job(form):
    from PROFIT_seq import env

    """Load form request from web UI"""
    job_name = form['jobName']
    start_time = int(form['startTime'])
    end_time = start_time + int(form['duration'])
    start_channel = int(form['startChannel'])
    end_channel = int(form['endChannel'])

    # Barcodes
    _barcode = form['barcode']
    if _barcode in ['classified', 'unclassified', 'all']:
        barcodes = _barcode
    else:
        barcodes = _barcode.split(',')

    # Target regions
    region_idx = defaultdict(list)
    for _target in form['target'].split('\n'):
        _target = _target.rstrip()
        if ":" not in _target:
            # Chr
            region_idx[_target] = form['mappedReads']
        else:
            if _target in env.Annotation:
                # Gene and transcript ID / Name
                chrom, st, en = env.Annotation[_target]
            else:
                # chr:start-end
                chrom, pos = _target.split(':')
                st, en = int(pos.split('-')[0]), int(pos.split('-')[1])
            if chrom in region_idx and isinstance(region_idx[chrom], str):
                continue
            region_idx[chrom].append((st, en, form['mappedReads']))
    
    region_idx["missed"] = form["missedReads"]
    region_idx["unmapped"] = form["unmappedReads"]
    region_idx["multi"] = form["multiMapped"]

    # Summarize job
    job = {
        'name': job_name,
        'time': (start_time, end_time),
        'ch': (start_channel, end_channel),
        'bc': barcodes,
        'index': region_idx,
        'group': form['mappedReads'],
    }
    return job


def load_gtf(gtf_file):
    gtf_idx = {}
    for parser in yield_gff(gtf_file):
        if parser.type == 'gene':
            gtf_idx[parser.gene_id] = [parser.contig, parser.start, parser.end]
            gtf_idx[parser.gene_name] = [parser.contig, parser.start, parser.end]
        elif parser.type == 'transcript':
            gtf_idx[parser.transcript_id] = [parser.contig, parser.start, parser.end]
            gtf_idx[parser.transcript_name] = [parser.contig, parser.start, parser.end]
        else:
            continue
    return gtf_idx


def hash_filename(fname):
    """Generate hash of file"""
    return hashlib.sha256(str(fname).encode()).hexdigest()


def to_bytes(bytes_or_str):
    """
    Return Instance of bytes
    """
    if isinstance(bytes_or_str, str):
        value = bytes_or_str.encode('utf-8')
    else:
        value = bytes_or_str
    return value


def to_str(bytes_or_str):
    """
    Return Instance of str
    """
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value


def grouper(iterable, n, fillvalue=None):
    """
    Collect data info fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    """

    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=None)


def flatten(x):
    """
    Flatten list of lists
    """
    import itertools

    flatted_list = list(itertools.chain(*x))
    return flatted_list


def check_dir(dir_name):
    """
    Create directory
    """
    dpath = Path(dir_name)
    if dpath.exists():
        if dpath.is_dir():
            pass
        else:
            sys.exit(f'Directory: {dpath}, conflict with existed files')
    else:
        dpath.mkdir()
    return dpath
