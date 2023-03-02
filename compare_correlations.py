#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import pyadjoint
import sys
from obspy import Trace
from obspy.geodetics import gps2dist_azimuth


def multitaper_traveltime(tr1, tr2, win, fqmin, fqmax):
    # get window in seconds
    win_t = win.copy()
    win_t[0] *= tr1.stats['delta']
    win_t[1] *= tr1.stats['delta']

    config = pyadjoint.Config(
        min_period=1.0/fqmax,
        max_period=1.0/fqmin,
        taper_percentage=0.1,
        measure_type='dt',
        min_cycle_in_window=3,
        dt_sigma_min=1.0,
        dlna_sigma_min=0.5,
        use_cc_error=True,
        use_mt_error=False,
        mt_nw=4.0,
        num_taper=5)

    adj_src = pyadjoint.calculate_adjoint_source(
        adj_src_type='cc_traveltime_misfit',
        config=config,
        observed=tr1,
        synthetic=tr2,
        window=[win_t],
        plot=False,
        adjoint_src=False)

    return adj_src.misfit


def compute_window(wf, lags, fqmin, fqmax):
    """
    Window surface waves on noise correlations

    :type syn: np array
    :type obs: np array
    """

    dt = abs(lags[1] - lags[0])
    nt = lags.size
    half_nt = ((nt - 1) / 2)

    # window center is the peak envelope
    # the envelope must be symmetric in time
    env = np.power(wf, 2)
    window_center = np.argmax(env)

    # window duration is <mult> times central period
    mult = 1.0
    T = (1./fqmax + 1./fqmin) / 2.
    window_offset = mult * T
    window_offset = int(window_offset / dt)

    # window index
    idx1 = max(window_center - window_offset, 0)
    idx2 = min(window_center + window_offset, nt - 1)

    # determine window for the positive branch
    if idx1 >= half_nt and idx2 > half_nt:
        window = [int(idx1), int(idx2)]
    elif idx1 < half_nt and idx2 <= half_nt:
        window = [int(2*half_nt - idx2), int(2*half_nt - idx1)]
    else:
        window = None

    return window


def compute_misft(waveforms1, lags1, waveforms2, lags2, fqmin, fqmax):
    nwaveforms = waveforms1.shape[0]
    windows = []
    t_misfit = []
    a_misfit = []

    for i in range(0, nwaveforms):
        wf1 = waveforms1[i, :]
        wf2 = waveforms2[i, :]

        # compute window
        win = compute_window(wf1, lags1, fqmin, fqmax)

        if win:
            # create obspy traces
            h1 = {'delta': abs(lags1[1] - lags1[0]), 'npts': lags1.size,
                  'channel': 'BXZ'}
            tr1 = Trace(wf1, header=h1)

            h2 = {'delta': abs(lags2[1] - lags2[0]), 'npts': lags2.size,
                  'channel': 'BXZ'}
            tr2 = Trace(wf2, header=h2)

            # compute traveltime misfit
            tmis = multitaper_traveltime(tr1, tr2, win, fqmin, fqmax)
            amis = 0.0
        else:
            tmis = np.nan
            amis = np.nan
            win = None

        windows.append(win)
        t_misfit.append(tmis)
        a_misfit.append(amis)

    return windows, t_misfit, a_misfit


def compute_az(rec, master_rec, rec_coor):
    net1, sta1 = rec.split('.')
    net2, sta2 = master_rec.split('.')

    coor1 = rec_coor[rec]
    coor2 = rec_coor[master_rec]

    dis, az, baz = gps2dist_azimuth(lat1=coor1[0],
                                    lon1=coor1[1],
                                    lat2=coor2[0],
                                    lon2=coor2[1])
    return dis, az


def read_stations_file(fname):
    dtype = [('sta', 'U20'), ('net', 'U20'), ('lat', 'f8'), ('lon', 'f8'),
             ('elv', 'f8'), ('z', 'f8')]

    X = np.loadtxt(fname, dtype=dtype)
    nrec = X.shape[0]
    rec_coor = {}

    for i in range(0, nrec):
        net = X['net'][i]
        sta = X['sta'][i]
        rec_coor[f'{net}.{sta}'] = [X['lat'][i], X['lon'][i]]

    return rec_coor


def read_ascii_waveform(fname):
    waveform = np.loadtxt(fname, dtype="float")
    time = waveform[:, 0]
    amp = waveform[:, 1]
    return time, amp


master_rec = 'XF.MA14'


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot one waveform.")
    parser.add_argument('path1')
    parser.add_argument('path2')
    parser.add_argument('file_ext')
    parser.add_argument('fqmin', type=float)
    parser.add_argument('fqmax', type=float)
    parser.add_argument('--stations', type=str, default=None)
    args = parser.parse_args()

    # scan data
    files = os.listdir(args.path1)
    files = [f for f in files if args.file_ext in f]
    nfiles = len(files)

    # get station-pair information
    if args.stations:
        rec_coor = read_stations_file(args.stations)

        pair_az = []
        pair_dis = []

        for f in files:
            net, sta, _, _ = f.split('.')
            rec = f'{net}.{sta}'
            dis, az = compute_az(rec, master_rec, rec_coor)
            pair_az.append(az)
            pair_dis.append(dis)

        pair_az = np.array(pair_az)
        pair_dis = np.array(pair_dis)

    # read data lags
    times1, _ = read_ascii_waveform(os.path.join(args.path1, files[0]))
    dt1 = abs(times1[1] - times1[0])
    ns1 = len(times1)
    maxlag = ((ns1 - 1) / 2) * dt1
    lags1 = np.linspace(-maxlag, maxlag, ns1)

    times2, _ = read_ascii_waveform(os.path.join(args.path2, files[0]))
    dt2 = abs(times2[1] - times2[0])
    ns2 = len(times2)
    maxlag = ((ns2 - 1) / 2) * dt2
    lags2 = np.linspace(-maxlag, maxlag, ns2)

    # read data
    waveforms1 = np.zeros((nfiles, ns1))
    waveforms2 = np.zeros((nfiles, ns2))

    for i, f in enumerate(files):
        _, a1 = read_ascii_waveform(os.path.join(args.path1, f))
        waveforms1[i, :] = a1

        _, a2 = read_ascii_waveform(os.path.join(args.path2, f))
        waveforms2[i, :] = a2

        # normalize
        waveforms2[i, :] = (a2 / np.max(np.abs(a2))) * np.max(np.abs(a1))

    # compute misfit
    if dt1 != dt2:
        print('Cannot compute misfit: waveforms have different sampling rate.')
        sys.exit()

    windows, t_misfit, a_misfit = compute_misft(waveforms1,
                                                lags1,
                                                waveforms2,
                                                lags2,
                                                fqmin=args.fqmin,
                                                fqmax=args.fqmax)

    # plot misfit
    t_misfit = np.array(t_misfit)

    plt.plot(t_misfit, 'o')
    plt.show()

    if args.stations:
        plt.plot(pair_dis, t_misfit, 'o')
        plt.show()

    print('total misfit ', np.nansum(t_misfit))

    # plot data
    offset = 0

    for i in range(0, waveforms1.shape[0]):
        if not windows[i]:
            print('no window ', files[i])
            continue

        w1 = waveforms1[i, :] + offset
        w2 = waveforms2[i, :] + offset
        win = windows[i]

        offset = max(w1.max(), w2.max())

        plt.plot(lags1, w1, color='k', alpha=0.8)
        plt.plot(lags2, w2, color='r', alpha=0.8)
        plt.plot(lags1[win[0]], w1[win[0]], 'ob')
        plt.plot(lags1[win[1]], w1[win[1]], 'ob')

    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.show()
