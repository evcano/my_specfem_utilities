#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram, windows
from obspy.signal.filter import bandpass

def read_ascii_waveform(file):
    waveform = np.loadtxt(file, dtype="float")
    time = waveform[:,0]
    amp = waveform[:,1]
    return time, amp


def plot_waveform(time, amp):
    plt.plot(time, amp, color="black")
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.show()


def plot_periodogram(amp, NS, DF):
    window = windows.get_window(window="hann", Nx=NS)

    freq, Pxx = periodogram(amp, fs=DF, window=window, detrend='constant',
                            return_onesided=True)
    plt.plot(freq, Pxx)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power')
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot one waveform.")
    parser.add_argument('file',help="Waveform file.")
    parser.add_argument('--f',type=float,nargs=2,default=None,help="Bandpass frequencies.")
    parser.add_argument('--p',action='store_true',help="Plot periodogram flag.")
    args = parser.parse_args()

    time, amp = read_ascii_waveform(args.file)
    DT = abs(time[1] - time[0])
    DF = 1.0 /  DT
    NS = len(time)
   
    if args.f:
       amp = bandpass(amp, freqmin=args.f[0], freqmax=args.f[1], df=DF,
                      corners=4, zerophase=True)

    if args.p == True:
        plot_periodogram(amp, NS, DF)
    else:
        plot_waveform(time, amp)
