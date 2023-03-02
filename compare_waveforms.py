#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pyadjoint
from scipy.signal import correlate, correlation_lags, periodogram, windows
from obspy.core import Trace


def Surf_Wave_MultiTaper(syn_tr, obs_tr, nt, dt, freqmin, freqmax, win):
    config = pyadjoint.Config(
        min_period=1./freqmax,
        max_period=1./freqmin,
        taper_percentage=0.0,
        measure_type='dt',
        min_cycle_in_window=3,
        dt_sigma_min=1.0,
        dlna_sigma_min=0.5,
        use_cc_error=True,
        use_mt_error=False,
        mt_nw=4.0,
        num_taper=5)

    adj_src = pyadjoint.calculate_adjoint_source(
        adj_src_type='multitaper_misfit',
        config=config,
        observed=obs_tr,
        synthetic=syn_tr,
        window=[win],
        plot=False,
        adjoint_src=True)

    return adj_src.misfit, adj_src.adjoint_source


def correlate_waveforms(tr1, tr2):
    corr = correlate(tr1.data, tr2.data)
    corr /= np.linalg.norm(tr1.data) * np.linalg.norm(tr2.data)
    corr_coef = np.max(corr)

    lags = correlation_lags(tr1.data.size, tr2.data.size, mode='full')
    max_lag = lags[np.argmax(corr)] * tr1.stats.delta
    return corr_coef, max_lag


def read_ascii_waveform(file):
    waveform = np.loadtxt(file, dtype="float")
    time = waveform[:, 0]
    amp = waveform[:, 1]
    return time, amp


def plot_waveform_overlap(time, amp1, amp2):
    difference = (amp1 - amp2)

    plt.subplot(211)
    plt.plot(time, amp1, color="black")
    plt.plot(time, amp2, color="red", linestyle='dashed')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')

    plt.subplot(212)
    plt.plot(time, difference, color='blue')
    plt.xlabel('Time (s)')
    plt.ylabel('Difference')

    plt.show()


def plot_periodogram_overlap(amp1, amp2, NS, DF):
    window = windows.get_window(window="hann", Nx=NS)

    freq1, Pxx1 = periodogram(amp1, fs=DF, window=window, detrend='constant',
                            return_onesided=True)
    freq2, Pxx2 = periodogram(amp2, fs=DF, window=window, detrend='constant',
                            return_onesided=True)
    difference = abs(Pxx1 - Pxx2)

    plt.subplot(211)
    plt.plot(freq1, Pxx1, color='black')
    plt.plot(freq2, Pxx2, color='red', linestyle='dashed')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power')

    plt.subplot(212)
    plt.plot(freq1, difference, color='blue')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Difference')

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot one waveform.")
    parser.add_argument('file1', help="Waveform file 1.")
    parser.add_argument('file2', help="waveform file 2.")
    parser.add_argument('--f', type=float, nargs=2, default=None,help="Bandpass frequencies.")
    parser.add_argument('--n', action='store_true', help="Normalize waveforms.")
    parser.add_argument('--p', action='store_true', help="Plot periodogram flag.")
    parser.add_argument('--r', action='store_true', help="Time reverse waveform 2.")
    args = parser.parse_args()

    time1, amp1 = read_ascii_waveform(args.file1)
    DT1 = abs(time1[1] - time1[0])
    DF1 = 1.0 / DT1
    NS1 = len(time1)

    time2, amp2 = read_ascii_waveform(args.file2)
    DT2 = abs(time2[1] - time2[0])
    DF2 = 1.0 / DT2
    NS2 = len(time2)

    DT1 = DT2
    if (DT1 != DT2):
        print(DT1, DT2)
        print("Both waveforms must have equal sampling rate.")
        sys.exit()

    if NS1 != NS2:
        dif = NS1 - NS2
        pad = np.abs(dif // 2)
        if dif > 0.0:
            amp2 = np.pad(amp2, (pad, pad), 'constant', constant_values=0.0)
            NS2 = NS1
        else:
            amp1 = np.pad(amp1, (pad, pad), 'constant', constant_values=0.0)
            NS1 = NS2

    tr1 = Trace(data=amp1, header={'delta': DT1, 'npts': NS1,
                                   'sampling_rate': DF1, 'channel': 'HXZ'})

    tr2 = Trace(data=amp2, header={'delta': DT2, 'npts': NS2,
                                   'sampling_rate': DF2, 'channel': 'HXZ'})


    if args.f:
        tr1.detrend('demean')
        tr1.detrend('linear')
        tr1.taper(0.1)
        tr1.filter('bandpass', freqmin=args.f[0], freqmax=args.f[1], corners=8,
                   zerophase=True)

        tr2.detrend('demean')
        tr2.detrend('linear')
        tr2.taper(0.1)
        tr2.filter('bandpass', freqmin=args.f[0], freqmax=args.f[1], corners=8,
                   zerophase=True)

    if args.n:
        tr1.normalize()
        tr2.normalize()

    if args.r:
        tr2.data = tr2.data[::-1]

    # print correlation coefficient and max lag
    corr_coef, max_lag = correlate_waveforms(tr1, tr2)
    print('Correlation coefficient: {}'.format(corr_coef))
    print('Maximum lag: {}'.format(max_lag))

    mute = False
    if mute:
        idx1 = int(26/0.05)
        idx2 = int(35/0.05)

        win = np.zeros(tr1.data.shape)
        win[idx1:idx2+1] += np.hanning(idx2+1-idx1)

        tr1.data *= win
        tr2.data *= win

    # plot figures
    if args.p:
        plot_periodogram_overlap(tr1.data, tr2.data,  NS1, DF1)

    plot_waveform_overlap(tr1.times(), tr1.data, tr2.data)

    #TODO: eliminate this
    dif = tr1.data - tr2.data
    dif = dif * 10.0E+20
    dif = dif[::-1]

    dif = dif.tolist()
    time2 = time2.tolist()

    out = np.array((time2, dif))
    np.savetxt('adjsrc', out.T)
