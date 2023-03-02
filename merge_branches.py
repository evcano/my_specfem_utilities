#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os


def read_ascii_waveform(filename):
    wf = np.loadtxt(filename, dtype='float')
    t = wf[:, 0]
    x = wf[:, 1]
    return t, x


def separate_branches(x):
    middle = (x.size + 1) // 2
    x_acaus = x[0:middle]
    x_caus = x[middle-1:]
    return x_acaus, x_caus


def merge_correlation(file1, file2):
    t1, c_ab = read_ascii_waveform(file1)
    t2, c_ba = read_ascii_waveform(file2)

    plt.plot(c_ab)
    plt.plot(c_ba)
    plt.show()

    if np.sum(t2 - t1) != 0.0:
        print('Error: The time axis of the correlations does not match.')
        return

    # time reverse due to symmetry C_ab(t) = C_ba(-t)
    c_ba = c_ba[::-1]


    c_ab_acaus, c_ab_caus = separate_branches(c_ab)
    c_ba_acaus, c_ba_caus = separate_branches(c_ba)

    # merge negative branch of c_ba with positive branch of c_ab
    corr = np.append(c_ba_acaus[:-1], c_ab_caus[:])

    return corr


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge two correlations')
    parser.add_argument('file1')
    parser.add_argument('file2')
    parser.add_argument('output1')
    parser.add_argument('output2')
    args = parser.parse_args()

    corr = merge_correlation(args.file1, args.file2)

    t, _ = read_ascii_waveform(args.file1)

    # save c_ab
    out = np.array([t, corr])
    out = out.T
    np.savetxt(args.output1, out, fmt='%1.8e')

    # save c_ba
    out[:, 1] = out[::-1, 1]
    np.savetxt(args.output2, out, fmt='%1.8e')
