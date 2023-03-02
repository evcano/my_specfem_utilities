#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import correlation_lags
from obspy.signal.filter import bandpass


def read_ascii_waveform(file):
    waveform = np.loadtxt(file, dtype="float")
    time = waveform[:, 0]
    amp = waveform[:, 1]
    return time, amp


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot one waveform.")
    parser.add_argument('datapath', help="Path to directory with data.")
    parser.add_argument('channel', help='Channel code.')
    parser.add_argument('extension', help='Data files extension.')
    parser.add_argument('--f', type=float, nargs=2, default=None, help="Bandpass frequencies.")
    parser.add_argument('--w', default=None, help="Path to waveform windows.")

    args = parser.parse_args()

    files = os.listdir(args.datapath)
    files = [f for f in files if args.channel in f]
    files = [f for f in files if args.extension in f]
    files.sort()

    offset = 0.

    for f in files:
        time, amp = read_ascii_waveform(os.path.join(args.datapath, f))

        net, sta, _, _ = f.split('.')
        stacode = '{}.{}'.format(net, sta)

        DT = abs(time[1] - time[0])
        DF = 1.0 / DT
        NS = len(time)

        NT = (NS + 1) / 2
        lags = correlation_lags(NT, NT, mode='full') * DT

        if args.f:
            amp = bandpass(amp, freqmin=args.f[0], freqmax=args.f[1], df=DF,
                           corners=4, zerophase=True)

        amp /= np.max(amp)
        amp += offset

        plt.plot(lags, amp, color='black', alpha=0.8)

        if args.w:
            win = np.load(os.path.join(args.w, f'{stacode}_window.npy'))
            plt.plot(lags[win[0]], offset, 'og')
            plt.plot(lags[win[1]], offset, 'og')

        offset = np.max(amp)

    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.show()
