#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
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
    parser.add_argument('--f', type=float, nargs=2, default=None,
                        help="Bandpass frequencies.")

    args = parser.parse_args()

    files = os.listdir(args.datapath)
    files = [f for f in files if args.channel in f]
    files = [f for f in files if args.extension in f]
    print('{} waveforms will be plotted'.format(len(files)))

    offset = 0.

    for f in files:
        wfpath = os.path.join(args.datapath, f)

        time, amp = read_ascii_waveform(wfpath)

        DT = abs(time[1] - time[0])
        DF = 1.0 / DT
        NS = len(time)

        if args.f:
            amp = bandpass(amp, freqmin=args.f[0], freqmax=args.f[1], df=DF,
                           corners=4)

        amp /= np.max(amp)
        amp += offset
        offset = np.max(amp)

        plt.plot(time, amp, color='black', alpha=0.8)

    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.show()
