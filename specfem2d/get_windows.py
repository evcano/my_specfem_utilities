import numpy as np
import os
from scipy.signal import correlation_lags


sta_dis = './reference_model/DATA/DISTANCE_TO_REFERENCE_STATION'

# data information
NT = 5000
DT = 0.1

# window parameters
vel = 3200.0  # m/s
dur = 40.0   # s

output_dir = './windows'

sta_dis = np.loadtxt(sta_dis, dtype=[('sta','U20'), ('dis','float')])
lags = correlation_lags(NT, NT, mode='full') * DT

windows = {}

for i, sta in enumerate(sta_dis['sta']):
    win_center = sta_dis['dis'][i] / vel
    win_start = max(win_center - dur / 2.0, 0.0)
    win_end = min(win_center + dur / 2.0, DT*NT)

    win_idx1 = np.argmin(np.abs(lags - win_start))
    win_idx2 = np.argmin(np.abs(lags - win_end))

    win = [win_idx1, win_idx2]
    np.save(os.path.join(output_dir, '{}_window.npy'.format(sta)), win)

print('finished computing windows')
