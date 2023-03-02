import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.fft as sf
from obspy.signal.filter import bandpass
from scipy.signal import correlation_lags, tukey


# PARAMETERS
# ==========
output_path = '.'

# spectrum information
spec_mean = 0.15
spec_std = 0.05

output_nt = 5000  # source function npts for one branch
output_dt = 0.1  # source function dt

# DONT EDIT BELOW THIS LINE
# =========================
freq = np.fft.rfftfreq(output_nt, d=output_dt)
spectrum = np.exp(-((freq - spec_mean)**2 / (2 * spec_std ** 2)))

taper = tukey(spectrum.size, 0.1)
spectrum *= taper

# get source time function
t = correlation_lags(output_nt, output_nt, mode='full') * output_dt
x = np.real(np.fft.irfft(spectrum, t.size, norm='backward'))
x = np.fft.fftshift(x)

# save source function
output = np.array([t, x], dtype='float32')
output = output.T
np.savetxt(os.path.join(output_path, 'S_squared'),
           output, fmt=['%1.7e', '%1.7e'])

# figures
plt.plot(freq, spectrum, 'k')
plt.show()

plt.plot(t, x)
plt.show()
