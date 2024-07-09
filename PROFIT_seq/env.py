import sys
import time
import mappy as mp
from traceback import format_exc
from pathlib import Path
from logging import getLogger
from threading import Thread
from collections import defaultdict
from read_until import ReadUntilClient
from pyguppy_client_lib.pyclient import PyGuppyClient

from PROFIT_seq.utils import Severity, send_message, load_gtf
Logger = getLogger("PROFIT_seq")

# Initialized server
Client = None
Caller = None
Aligner = None
Debugging = False
Annotation = {}
RedFlag = False
Sequence = {}
Jobs = []

# Init / Waiting / Running / Stopped / Finished
Worker = {"status": "Init", "worker": None, "message": "Initializing"}
ProtocolRun = {}

# Storage reads for each barcode / target
Barcode_coverage = defaultdict(int)
Target_coverage = defaultdict(int)
Target_blacklist = {}
Unblock_header = ['read_id', 'channel', 'read_number', 'seq_length', 'action', 'reason', 'alignment', 'sequence']
Log_header = ['time', 'total', 'get_read', 'basecalling', 'alignment', 'judge', 'action', 'update']

def initializer(minknow_host, minknow_port, guppy_address, guppy_config, mm_idx, gtf_file):
    """Initialize servers"""
    global Caller, Client, Aligner, Annotation, ProtocolRun, Sequence, Debugging
    Logger.debug('Initializing PROFIT-seq runner ...')

    # Step1. Connect to MinKNOW server
    Logger.debug(f"MinKNOW server: {minknow_host}:{minknow_port}")
    try:
        read_until_client = ReadUntilClient(
            mk_host=minknow_host, mk_port=minknow_port, one_chunk=False,
            # cache_type=AccumulatingCache,
            # filter_strands=True,
            # calibrated_signal=False
        )
        batch_size = read_until_client.channel_count
        Logger.debug(f"Using batch size of {batch_size}")
    except Exception as e:
        Logger.error("Failed to connect to MinKNOW server. Please make sure the MinION device and sequencing Flow Cell is properly installed!")
        Logger.error(format_exc())
        Logger.error(sys.exc_info()[2])
        sys.exit(1)

    # assert read_until_client.connection.acquisition.current_status().status == MinknowStatus.READY
    send_message(
        read_until_client.connection,
        "Message from FullSpeed: Successfully connected!",
        Severity.INFO,
    )
    Client = read_until_client
    ProtocolRun['device_id'] = Client.connection.device.get_device_info().device_id
    Logger.info("Connected to MinKNOW server")

    # Step2. Connect to ont-guppy basecaller
    Logger.debug(f"Guppy server: {guppy_address}, using config: {guppy_config}")
    caller = PyGuppyClient(
        address=guppy_address,
        config=guppy_config,
        server_file_load_timeout=180,
    )
    caller.connect()
    caller.disconnect()
    Caller = {"address": guppy_address, "config": guppy_config, "server_file_load_timeout": 180}
    Logger.info("Connected to Guppy server")

    # Step3. Connect to minimap2 aligner
    Logger.debug(f"Loading reference index {mm_idx}")
    # aligner = None
    aligner = mp.Aligner(mm_idx)
    Aligner = aligner
    Logger.info("Loaded reference index")

    if gtf_file is not None:
        annotation = load_gtf(gtf_file)
        Annotation = annotation
        Logger.info("Loaded annotation gtf")

    Worker['status'] = 'Ready'
    Worker['message'] = 'Ready for sequencing'

    return 0


def split_segments(segments, is_closed=False):
    """Split segments into sub intervals"""
    if len(segments) == 0:
        return [], {}

    if len(segments) == 1:
        return segments, {segments[0]: [0, ]}

    nodes = []
    for i, j in segments:
        if is_closed:
            nodes += [i-1, j]
        else:
            nodes += [i, j]
    nodes = sorted(list(set(nodes)))
    edges = [[i, j] for i, j in zip(nodes[:-1], nodes[1:])]

    segments_d = defaultdict(list)
    for x in sorted(list(set(segments))):
        p = x[0]-1 if is_closed else x[0]
        q = x[1]
        for i, (l_st, l_en) in enumerate(edges):
            if q <= l_st or l_en <= p:
                continue
            segments_d[x].append(i)

    return edges, segments_d


def get_jobs():
    time_sections = [job['time'] for job in Jobs]
    time_nodes, time_d = split_segments(time_sections, False)

    ch_sections = [job['ch'] for job in Jobs]
    ch_nodes, ch_d = split_segments(ch_sections, True)

    job_mtx = defaultdict(dict)
    for job in Jobs:
        for i in time_d[job['time']]:
            for j in ch_d[job['ch']]:
                job_mtx[i].setdefault(j, []).append(job)

    return time_nodes, ch_nodes, job_mtx

