import os
import time
import traceback
from pathlib import Path
from logging import getLogger
from threading import Thread
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pyguppy_client_lib.pyclient import PyGuppyClient
from minknow_api.acquisition_pb2 import MinknowStatus

from PROFIT_seq import env
from PROFIT_seq.env import get_jobs
from PROFIT_seq.utils import Severity, send_message, basecall
Logger = getLogger("PROFIT_seq")


def worker_callback(fn):
    if fn.cancelled():
        Logger.error(f"{Path(__file__).name} - worker cancelled")
    elif fn.done():
        error = fn.exception()
        if error:
            Logger.error(f"{Path(__file__).name} - {type(fn.exception()).__name__} in worker thread: {error.__cause__}")
            Logger.error(traceback.format_exc())
        else:
            Logger.info(fn.result())


def start_unblock():
    try:
        Logger.info("Running unblock jobs ...")

        # Merge unblock jobs
        time_nodes, ch_nodes, job_mtx = get_jobs()

        env.Worker['status'] = 'Waiting'
        env.Worker['message'] = 'Waiting for connection'

        # Wait for job start
        while env.Client.connection.acquisition.current_status().status != MinknowStatus.PROCESSING:
            Logger.warning('No processing protocol detected, waiting for retry in 10 seconds ...')
            time.sleep(10)
            if env.RedFlag:
                return 0
        unblock_start = time.time()

        env.Client.run()
        send_message(
            env.Client.connection,
            "Message from FullSpeed: Start unblocking jobs!",
            Severity.INFO,
        )

        # Load run info
        run_info = env.Client.connection.protocol.get_current_protocol_run()
        param = {i.split('=')[0].strip('-'): i.split('=')[1] for i in run_info.args if '=' in i}
        env.ProtocolRun.update({
            'device_id': run_info.device.device_id,
            'product': run_info.flow_cell.user_specified_product_code,
            'flow_cell': run_info.flow_cell.user_specified_flow_cell_id,
            'experiment': run_info.user_info.protocol_group_id.value,
            'sample': run_info.user_info.sample_id.value,
            'start_time': run_info.start_time.seconds,
            'unblock_start': unblock_start,
            'output_path': run_info.output_path,
            'barcode_kits': run_info.meta_info.tags['barcoding kits'],
            'run_dir': run_info.output_path,
            'lock_file': Path(run_info.output_path) / 'fusion.lock',
            'unblock_file': Path(run_info.output_path) / 'unblock_reads.txt',
            'log_file': Path(run_info.output_path) / 'unblock_log.txt',
            'channel_count': run_info.flow_cell.channel_count,
            'experiment_time': param['experiment_time'],
            'barcoding_kits': [i.strip("'\"[]") for i in param['barcoding_kits'].split(',')] if 'barcoding_kits' in param else [],
        })
        env.ProtocolRun['lock_file'].write_text('0')
        env.ProtocolRun['unblock_fh'] = open(env.ProtocolRun['unblock_file'], 'a')
        env.ProtocolRun['log_fh'] = open(env.ProtocolRun['log_file'], 'a')
        if env.ProtocolRun['unblock_file'].stat().st_size == 0:
            env.ProtocolRun['unblock_fh'].write('\t'.join(env.Unblock_header) + '\n')
        if env.ProtocolRun['log_file'].stat().st_size == 0:
            env.ProtocolRun['log_fh'].write('\t'.join(env.Log_header) + '\n')

        Logger.info(
            f'Device: %s, Flow cell (%s): %s, Experiment: %s, Sample: %s',
            env.ProtocolRun['device_id'], env.ProtocolRun['product'], env.ProtocolRun['flow_cell'],
            env.ProtocolRun['experiment'], env.ProtocolRun['sample'],
        )

        # Connect to guppy caller
        if 'barcoding_kits' in param:
            caller = PyGuppyClient(
                address=env.Caller['address'],
                config=env.Caller['config'],
                barcode_kits=' '.join(env.ProtocolRun['barcoding_kits']),
                server_file_load_timeout=env.Caller['server_file_load_timeout'],
            )
        else:
            caller = PyGuppyClient(
                address=env.Caller['address'],
                config=env.Caller['config'],
                server_file_load_timeout=env.Caller['server_file_load_timeout'],
            )
        caller.connect()
        env.Caller = caller

        with ThreadPoolExecutor(max_workers=1) as executor:
            break_time = 0.4
            worker = executor.submit(run_unblock, time_nodes, ch_nodes, job_mtx, break_time, 10, break_time)
            worker.add_done_callback(worker_callback)

        env.ProtocolRun['lock_file'].write_text('1')
        env.ProtocolRun['lock_file'].unlink()
        env.ProtocolRun['unblock_fh'].close()
        env.ProtocolRun['log_fh'].close()
        env.Client.reset()
        env.Caller.disconnect()
        Logger.info("Finished")

    except KeyboardInterrupt:
        send_message(
            env.Client,
            'Message from FullSpeed: Stopped adaptive sequencing due to keyboard interupt.',
            Severity.WARN,
        )
    return 0


def run_unblock(time_nodes, ch_nodes, job_mtx, throttle=0.4, retry_interval=10, unblock_duration=0.4):
    env.ProtocolRun['lock_file'].write_text(str(os.getpid()))
    env.Worker['status'] = 'Running'
    env.Worker['message'] = 'Running'

    due_time = env.ProtocolRun['start_time'] + float(env.ProtocolRun['experiment_time']) * 60 * 60
    while env.Client.is_running and env.ProtocolRun['lock_file'].exists() and env.Worker['status'] != 'Stopped':
        t0 = time.time()
        if t0 > due_time:
            break
        if env.RedFlag:
            return 0

        # Get batch reads
        read_batch = env.Client.get_read_chunks(batch_size=env.ProtocolRun['channel_count'], last=True)
        t1 = time.time()

        # Basecall reads
        try:
            called_batch = basecall(
                guppy_client=env.Caller,
                reads=read_batch,
                dtype=env.Client.signal_dtype,
                daq_values=env.Client.calibration_values
            )
            called_batch = list(called_batch)
        except ConnectionError as e:
            Logger.warning(f'Failed to perform reads basecalling, skip this batch ...\n{e}')
            batch_time = time.time() - t0
            if batch_time < throttle:
                time.sleep(throttle - batch_time)
            continue

        t2 = time.time()

        # Get current time
        tmp_blacklist = []
        tmp_whitelist = []
        tmp_waitlist = []
        batch_lines = []
        n_black, n_white, n_wait = 0, 0, 0

        ctime = (t0 - env.ProtocolRun['start_time']) / 60
        _t = [x for x, (p, q) in enumerate(time_nodes) if p < ctime <= q]
        t_align = 0
        for (channel, read_number), read in called_batch:
            # Filter using time and channel number
            _ch = [x for x, (p, q) in enumerate(ch_nodes) if p < channel <= q]
            if len(_t) == 0 or len(_ch) == 0:
                tmp_whitelist.append((channel, read_number))
                continue

            assert len(_t) == 1 and len(_ch) == 1

            read_id = read[0]['metadata']['read_id']
            read_bc = read[0]['metadata']['barcode_arrangement']
            read_seq = read[0]['datasets']['sequence']
            if read_id in env.Sequence:
                read_seq = env.Sequence[read_id] + read_seq

            action, reason, hits, tmp_t = anubis(read_id, read_seq, read_bc, job_mtx[_t[0]][_ch[0]])
            t_align += tmp_t

            if action in ['unblock', 'balance_unblock']:
                if read_id in env.Sequence:
                    del env.Sequence[read_id]
                tmp_blacklist.append((channel, read_number))
                n_black += 1
            elif action in ['stop_receiving', 'balance_stop_receiving']:
                if read_id in env.Sequence:
                    del env.Sequence[read_id]
                tmp_whitelist.append((channel, read_number))
                n_white += 1
            elif action == "wait":
                # Store tmp sequences
                env.Sequence[read_id] = read_seq
                tmp_waitlist.append((channel, read_number))
                n_wait += 1

            tmp_hit = None if hits is None else ",".join([f"{x.ctg}:{x.r_st}-{x.r_en}" for x in hits])
            tmp_line = [
                read_id,
                channel,
                read_number,
                len(read_seq),
                action,
                reason,
                tmp_hit,
            ]
            batch_lines.append(tmp_line)

        t3 = time.time()
        if tmp_blacklist:
            env.Client.unblock_read_batch(tmp_blacklist, duration=unblock_duration)
            env.Client.stop_receiving_batch(tmp_blacklist)

        if tmp_whitelist:
            env.Client.stop_receiving_batch(tmp_whitelist)

        if batch_lines:
            for tmp_line in batch_lines:
                env.ProtocolRun['unblock_fh'].write('\t'.join([str(x) for x in tmp_line]) + '\n')
        t4 = time.time()

        # Wait till throttle
        batch_time = time.time() - t0
        Logger.info(f"Unblock {n_black}/Passed {n_white}/Proceed {n_wait} in {batch_time}s, {len(env.Sequence)} cached")

        # Update coverage for balance jobs
        update_coverage()
        t5 = time.time()

        env.ProtocolRun['log_fh'].write(
            '\t'.join([str(x) for x in [t5, t5-t0, t1-t0, t2-t1, t_align, t3-t2-t_align, t4-t3, t5-t4]]) + '\n'
        )

        batch_time = time.time() - t0
        if batch_time < throttle:
            time.sleep(throttle - batch_time)

        # Retry if lost connection
        retry_time = 0
        while env.Client.connection.acquisition.current_status().status != MinknowStatus.PROCESSING:
            retry_time += 1
            Logger.warning('Lost connection, retrying for the %s time ...', retry_time)
            time.sleep(retry_interval)
            if retry_time > 10:
                Logger.info(f'No connection for {retry_time*retry_interval:d} seconds, abort')
                break
        else:
            continue
        break

    env.Worker['status'] = 'Finished'
    env.Worker['message'] = 'Unblock finished'
    return 0


def anubis(read_id, read_seq, read_bc, jobs):
    """Judge whether continue sequencing for specific read"""
    bc_jobs = []
    for job in jobs:
        if job['bc'] == 'all':
            bc_jobs.append(job)
        elif job['bc'] == 'classified':
            if read_bc != 'unclassified':
                bc_jobs.append(job)
        elif job['bc'] == 'unclassified':
            if read_bc == 'unclassified':
                bc_jobs.append(job)
        else:
            if read_bc in job['bc']:
                bc_jobs.append(job)

    # No applicable filters
    if len(bc_jobs) > 1:
        Logger.error(f"Conflict jobs: {bc_jobs}")
        return "stop_receiving", "barcode_conflict", None, 0
    if len(bc_jobs) == 0:
        return "stop_receiving", "barcode_nojob", None, 0

    job = bc_jobs[0]
    if 'all' in job['index']:
        return job['index']['all'], "barcode_all", None, 0

    # Minimap2 Alignment filter
    t0 = time.time()
    hits = [hit for hit in env.Aligner.map(read_seq) if hit.is_primary]
    t1 = time.time()

    # contigs = set([i.ctg for i in hits])

    # Unmapped
    if len(hits) == 0:
        # Disable when testing timing
        # return job['index']['unmapped'], "unmapped", None, t1 - t0
        if len(read_seq) > 500: # TODO: determine minimum unmapped length
            return job['index']['unmapped'], "unmapped", None, t1-t0
        else:
            return "wait", "unmapped_too_short", None, t1-t0

    # Mapped
    if 'mapped' in job['index']:
        return job['index']['mapped'], "mapped", hits, t1-t0

    # Multi-mapped
    if len(hits) > 1:
        if 'multi' in job['index']:
            return job['index']['multi'], "multi", hits, t1-t0

    actions_ = defaultdict(list)
    for hit in hits:
        contig, st, en = hit.ctg, hit.r_st, hit.r_en
        if contig in job['index']:
            for r_st, r_en, action in job['index'][contig]:
                if r_en < st:
                    continue
                elif en < r_st:
                    break
                else:
                    actions_[action].append(f"{contig}:{r_st}-{r_en}")

    if len(actions_) == 0:
        return job['index']['miss'], "miss", hits, t1-t0

    if 'unblock' in actions_:
        return 'unblock', ",".join(actions_['unblock']), hits, t1-t0

    if 'stop_receiving' in actions_:
        return 'stop_receiving', ",".join(actions_['stop_receiving']), hits, t1-t0

    if 'wait' in actions_:
        return 'wait', ",".join(actions_['wait']), hits, t1-t0

    if 'balance' in actions_:
        balance_actions_ = []
        for x in actions_['balance']:
            if env.Target_blacklist and x in env.Target_blacklist:
                balance_actions_.append("unblock")
            else:
                env.Target_coverage[x] += 1
                balance_actions_.append("stop_receiving")
        if 'unblock' in balance_actions_:
            return "balance_unblock", ",".join(actions_['balance']), hits, t1-t0
        if "stop_receiving" in balance_actions_:
            return "balance_stop_receiving", ",".join(actions_['balance']), hits, t1-t0

    Logger.error(f"Unkown operation: {actions_}")

    return "stop_receiving", ",".join(actions_[list(actions_)[0]]), hits, t1-t0


def update_coverage():
    env.Target_blacklist = {}
    sorted_keys = sorted(env.Target_coverage, key=lambda x: env.Target_coverage[x], reverse=True) # From inf->0
    if len(sorted_keys) > 0:
        idx = int(len(sorted_keys) * .2)
        idx = max(idx, 1)
        for i in sorted_keys[:idx]:
            env.Target_blacklist[i] = 1
    return 0
