import os
import time
import json
import socket
from logging import getLogger
from threading import Thread
from flask import Flask, render_template, url_for, request, jsonify
from PROFITseq import env
from PROFITseq.unblock import start_unblock
from PROFITseq.utils import load_toml, load_prim_job
from app.plot import plot_jobs

Logger = getLogger("PROFITseq")
app = Flask(__name__, template_folder='templates', static_folder='static')


@app.route('/')
def index():
    hostname = socket.gethostname()
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(('8.8.8.8', 80))
        local_ip = s.getsockname()[0]
    except:
        local_ip = socket.gethostbyname(hostname + ".local")
    finally:
        if s:
            s.close()
    return render_template(
        'index.html', hostname=hostname, local_ip=local_ip, device=env.ProtocolRun['device_id']
    )


@app.route('/plot', methods=['POST'])
def plot():
    graphJSON = plot_jobs()
    return jsonify({"runStatus": env.Worker['status'], "graph": graphJSON})


@app.route('/submit',  methods=['GET', 'POST'])
def submit():
    Logger.info(request.form)
    Logger.info(request.files)

    if request.files['jobFile'].filename != '':
        jobs = load_toml(request.files['jobFile'])
        env.Jobs += jobs
    else:
        job = load_prim_job(request.form)
        env.Jobs.append(job)

    return jsonify({"status": "Success"})


@app.route('/clear')
def clear():
    """Clear existing jobs"""
    env.Jobs = []
    return jsonify({"status": "Cleared"})


@app.route('/run')
def run():
    """Start unblock process when query"""
    thread = Thread(target=start_unblock)
    thread.start()
    env.Worker['worker'] = thread
    return jsonify({"status": "Running"})


@app.route('/stop')
def stop():
    """Stop unblock process when query"""
    env.RedFlag = True
    env.Worker['message'] = "Waiting for experiment stop"
    if env.Worker['worker'] is not None:
        env.Worker['worker'].join()
    env.Worker['status'] = 'Stopped'
    env.Worker['message'] = 'Unblock stopped'
    return jsonify({"status": "Stopped"})


@app.route('/prog')
def progress():
    """Update sequencing progress when query"""
    from PROFITseq.utils import convert_time
    ptext = env.Worker['message']
    dtext = env.ProtocolRun['device_id']
    stext = env.ProtocolRun['sample'] if 'sample' in env.ProtocolRun else "TBD"
    pnum = 0

    # Init / Ready (Connected to MinKNOW) / Running (Running Unblock Jobs) / Finished (Job finished)
    if env.Worker['status'] in ['Ready', 'Waiting']:
        pass
    elif env.Worker['status'] == 'Running':
        duration = time.time() - float(env.ProtocolRun['start_time'])
        experiment_time = float(env.ProtocolRun['experiment_time'])
        ptext = f"{convert_time(duration)} / {experiment_time}h"
        pnum = duration / 60 / 60 / experiment_time * 100
    elif env.Worker['status'] in ['Finished', 'Stopped']:
        pnum = 100
    else:
        Logger.info(f"Unknown status: {env.Worker['status']}")

    return jsonify({"runStatus": env.Worker['status'], "progNumber": pnum, "progText": ptext, "deviceText": dtext, "sampleText": stext})
