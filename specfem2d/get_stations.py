import matplotlib.pyplot as plt
import numpy as np


net_code = 'DB'
sta_code = ['A', 'B', 'C', 'D', 'E']

# reference station coordinates
x0 = 0.0
y0 = 0.0

# depth and burial coordinates of all stations
z = 0.0
b = 0.0

# generate stations coordinates
min_wavelength = 6.6 * 3200.0
mult = [2.0, 3.0, 4.0, 5.0, 6.0]

angle = [270.0, 300.0, 330.0, 0.0, 30.0, 60.0, 90.0]
angle = np.deg2rad(angle)


# write file
_f = open('STATIONS', 'w')
_f2 = open('DISTANCE_TO_REFERENCE_STATION', 'w')

j = 0
i = 0

sta = '{}{:03d}'.format(sta_code[j], i)

_f.write('{} {} {:.4f} {:.4f} {:.4f} {:.4f}\n'.format(
    sta, net_code, y0, x0, z, b))

_f2.write('{}.{} {}\n'.format(net_code, sta, 0.0))

for j, r in enumerate(mult):
    for i, a in enumerate(angle):
        i += 1

        dis = r * min_wavelength
        x = (dis * np.cos(a)) + x0
        y = (dis * np.sin(a)) + y0
        plt.plot(x, y, 'o')

        sta = '{}{:03d}'.format(sta_code[j], i)

        _f.write('{} {} {:.4f} {:.4f} {:.4f} {:.4f}\n'.format(
            sta, net_code, x, y, z, b))

        _f2.write('{}.{} {}\n'.format(net_code, sta, dis))

        theta = (np.arctan2(y-y0, x-x0) + 2*np.pi) % (2*np.pi)

_f.close()
_f2.close()

plt.plot(x0, y0, 'o')
plt.show()
