import time
import json
import numpy as np
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
from flask import Flask, render_template

from PROFITseq import env

Action_color = {
    "unblock": "#a53f97", "stop_receiving": "#f57e20", "balance": "#90c6e7", "unblock all": "#17bb75",
}


def plot_jobs():
    xmax = 72
    xprog = 0
    ymax = 512

    if env.Jobs:
        xmax = np.max([job['time'][1]/60 for job in env.Jobs])

    # Get run info
    if 'experiment_time' in env.ProtocolRun:
        xmax = float(env.ProtocolRun['experiment_time'])
        xprog = (time.time() - env.ProtocolRun['start_time'])/60/60

    if xmax >= 72:
        xdtick = 12
    elif xmax >= 48:
        xdtick = 6
    elif xmax >= 24:
        xdtick = 4
    elif xmax >= 12:
        xdtick = 3
    elif xmax >= 6:
        xdtick = 1
    elif xmax >= 3:
        xdtick = 0.5
    else:
        xdtick = int(10 * xmax / 10) / 10

    # Consruct plot
    fig = go.Figure()

    if env.Jobs:
        n_actions = len(set([x['group'] for x in env.Jobs]))
        for job in env.Jobs:
            x1, x2 = job['time'][0]/60, job['time'][1]/60
            y1, y2 = job['ch']
            c = Action_color[job['group']]

            # t = f"barcode: {job['bc']}"
            # if 'target' in job:
            #     t += f", target: {job['target']}"
            # t += f", action: {job}"
            t = job['name']

            if n_actions > 1:
                fig.add_trace(go.Scatter(
                    x=[x1,x1,x2,x2,x1], y=[y1,y2,y2,y1,y1],
                    mode="none", fill='toself', fillcolor=c, opacity=0.4,
                    name=job['name'], text=job['name'], hoverinfo="none", showlegend=False,
                ))
            else:
                fig.add_trace(go.Scatter(
                    x=[x1,x1,x2,x2,x1], y=[y1,y2,y2,y1,y1],
                    mode="none", fill='toself', opacity=0.4,
                    name=job['name'], text=job['name'], hoverinfo="none", showlegend=False,
                ))
            fig.add_trace(go.Scatter(
                x=[(x1+x2)/2, ], y=[(y1+y2)/2, ],
                opacity=0.4,
                mode='markers', marker=dict(size=[0, ]),
                hoverinfo='text', text=[t,], showlegend=False,
            ))

    # Sequencing progress
    fig.add_vrect(
      x0=0, x1=xprog, fillcolor="#bdbdbd", opacity=0.5, layer="above", line_width=0,
    )

    fig.update_xaxes(range=[0, xmax], fixedrange=True)
    fig.update_layout(
        hovermode="x",
        hoverlabel=dict(
            bgcolor="white",
            font_size=14,
            font_family="Nunito",
        ),
        plot_bgcolor='#f6f7fa',
        margin=go.layout.Margin(l=0, r=0, b=0, t=0),
        xaxis = dict(
            fixedrange = True, range = [0, xmax],
            tickmode = 'linear', tick0 = 0, dtick = xdtick,
            title = 'Time (h)', title_font = dict(size=14, family='Nunito'),
        ),
        yaxis = dict(
            fixedrange = True, range = [0, ymax],
            tickmode = 'array', tickvals = [0, 128, 256, 384, 512], ticktext = ["0 ", "128 ", "256 ", "384 ", "512 "],
            title = "Channel Number", title_font = dict(size=14, family='Nunito'),
        )
    )

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return graphJSON
