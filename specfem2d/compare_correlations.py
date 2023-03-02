#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import correlation_lags
from matplotlib.patches import Rectangle
from obspy.core import Trace
from obspy.signal.filter import bandpass


def compute_misfit(time, amp1, amp2, window, misfit_type, fqband):
    DT = abs(time[1] - time[0])
    DF = 1.0 / DT
    NT = time.size

    header = {'delta': DT, 'npts': NT, 'sampling_rate': DF, 'channel': 'BXY'}
    tr1 = Trace(amp1, header=header)
    tr2 = Trace(amp1, header=header)

    misfit = [np.nan, np.nan]

    for i in range(0, 2):
        # causal branch
        if i == 0:
            window = window
        # acausal branch
        elif i == 1:
            window = [NT - 1 - win[1], NT - 1 - win[0]]

        conf = pyadjoint.Config(
            min_period=1.0/fqband[1],
            max_period=1.0/fqband[0],
            taper_percentage=0.05,
            measure_type='dt',
            min_cycle_in_window=3,
            dt_sigma_min=1.0,
            dlna_sigma_min=0.5,
            use_cc_error=True,
            use_mt_error=False,
            mt_nw=4.0,
            num_taper=5)

        adjsrc = pyadjoint.calculate_adjoint_source(
            adj_src_type=misfit_type,
            config=conf,
            observed=tr1,
            synthetic=tr2,
            window=window,
            plot=False,
            adjoint_src=False)

        misfit[i] = adjsrc.misfit

    return misfit


def read_ascii_waveform(file):
    waveform = np.loadtxt(file, dtype="float")
    time = waveform[:, 0]
    amp = waveform[:, 1]
    return time, amp


data_path1 = './reference_model/OUTPUT_STEP_2'
data_path2 = './coast_model/OUTPUT_STEP_2'
channel = 'BXY'
extension = 'semd'
fqband = None
windows_path = './windows'
compute_misfit = ''
misfit_file = None

# figure related
figname = None
plt_maxlag = 100.0
plt_fsize = 16.0
plt_xinch = 10.0
plt_yinch = 10.0

if __name__ == "__main__":
    files = os.listdir(data_path1)
    files = [f for f in files if f'{channel}.{extension}' in f]
    files.sort()

    fig, ax = plt.subplots()
    fig.set_size_inches(plt_xinch, plt_yinch)
    offset = 0.

    sta_misfit = {}

    for f in files:
        # read data
        net, sta, _, _ = f.split('.')
        stacode = '{}.{}'.format(net, sta)

        time, amp1 = read_ascii_waveform(os.path.join(data_path1, f))
        _, amp2 = read_ascii_waveform(os.path.join(data_path1, f))

        DT = abs(time[1] - time[0])
        DF = 1.0 / DT
        NT = len(time)
        BRANCH_NT = (NT - 1) // 2

        lags = correlation_lags(BRANCH_NT, BRANCH_NT, mode='full') * DT

        if windows_path:
            win = np.load(os.path.join(windows_path, f'{stacode}_window.npy'))

        # filter data
        if fqband:
            amp1 = bandpass(amp1, freqmin=fqband[0], freqmax=fqband[1], df=DF,
                            corners=4, zerophase=True)

            amp2 = bandpass(amp2, freqmin=fqband[0], freqmax=fqband[1], df=DF,
                            corners=4, zerophase=True)

        # normalize data
        amp1 /= np.max(amp1)
        amp2 /= np.max(amp2)

        # compute misfit
        if windows_path and compute_misfit:
            misfit = compute_misfit(time, amp1, amp2, win, compute_misfit)
            sta_misfit[stacode] = misfit.copy()

        # ==============================================================================
        # plotting
        # ==============================================================================
        amp1 += offset
        amp2 += offset

        ax.plot(lags, amp1, color='black', alpha=0.8)
        ax.plot(lags, amp2, color='blue', alpha=0.8)

        if windows_path:
            # causal branch
            ax.add_patch(Rectangle((lags[win[0]], amp1.min()),
                                   lags[win[1]]-lags[win[0]],
                                   amp1.max()-amp1.min(), color='k', alpha=0.1))
            # acausal branch
            win = [2*BRANCH_NT-win[1], 2*BRANCH_NT-win[0]]
            ax.add_patch(Rectangle((lags[win[0]], amp1.min()),
                                   lags[win[1]]-lags[win[0]],
                                   amp1.max()-amp1.min(), color='k', alpha=0.1))

        offset = np.max(amp1)

    plt.xlabel('Time (s)', fontsize=plt_fsize)
    plt.ylabel('Amplitude', fontsize=plt_fsize)
    plt.xticks(fontsize=plt_fsize*0.8)
    plt.yticks(fontsize=plt_fsize*0.8)
    plt.xlim((-plt_maxlag, plt_maxlag))
    plt.grid()

    if not figname:
        plt.show()
    else:
        print('saving figure')
        plt.savefig(figname, dpi=600, bbox_inches='tight', facecolor='white')

    if misfit_file:
        print('saving misfit')
