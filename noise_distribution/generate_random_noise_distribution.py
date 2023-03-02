#!/usr/bin/env python3
#
# Generates noise distribution and direction binary files
# according to parfile_noise.yaml.
#
# The binary files contain float values for all GLL points
# on the mesh free-surface.
#
# Make sure the free-surface corresponds to the mesh surface
# where you want noise sources to be located.
#
# IMPORTANT: This script needs proc*_free_surface.vtk files in
# DATABASES_MPI. To generate them, set SAVE_MESH_FILES = .true.
# in DATA/Par_file and run xgenerate_databases.
#


import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import yaml
from convert_utm2geo import utm_geo
from scipy.interpolate import RegularGridInterpolator


def fft_idx_gen(size):
    '''
    generate fft indexes
    '''

    a = np.arange(0, size//2+1)
    b = np.arange(1, size//2)
    b = b[::-1] * -1
    kidx = np.append(a, b)
    return kidx


def gaussian_spectrum(kx, ky, std):
    mean = 0.0
    k = np.sqrt(kx**2 + ky**2)
    amp = np.exp(-(k-mean)**2/(2*std**2))
    return amp


def gaussian_random_field(nx, dx, wavelength):
    std = 1.0 / (2.0*wavelength)

    kres = 1.0 / (nx * dx)
    kidx = fft_idx_gen(nx)
    kaxis = kidx * kres

    noise = np.random.normal(size=(nx, nx))
    noise_fft = np.fft.fftn(noise)

    amplitude = np.zeros((nx, nx))
    for i, kx in enumerate(kaxis):
        for j, ky in enumerate(kaxis):
            amplitude[i, j] = gaussian_spectrum(kx, ky, std)

    random_field = np.real(np.fft.ifftn(noise_fft * amplitude))
    return random_field


def gaussian_distribution(distribution, mask_noise, xcoor, ycoor):
    if distribution['utm_zone']:
        center_x, center_y = utm_geo(distribution['center_x'],
                                     distribution['center_y'],
                                     distribution['utm_zone'],
                                     2)
    else:
        center_x = distribution['center_x']
        center_y = distribution['center_y']

    dist = np.zeros((xcoor.size))
    a = np.array([center_x, center_y])

    for i in range(0, dist.size):
        b = np.array([xcoor[i], ycoor[i]])
        dist[i] = np.linalg.norm(a-b)

    gaussian = np.exp(- (dist ** 2) / (2 * distribution['sigma_m'] ** 2))
    mask_noise += gaussian * distribution['weight']
    return mask_noise


def uniform_distribution(distribution, mask_noise):
    mask_noise[:] += distribution['weight']
    return mask_noise


def ocean_distribution(distribution, mask_noise, zcoor):
    ocean_idx = np.where(zcoor < 0.0)[0]

    if not ocean_idx.any():
        print('Skipping ocean distribution: There are no gll points on\
               the ocean.')
    else:
        mask_noise[ocean_idx] += distribution['weight']
    return mask_noise


def read_gll_coordinates(ifile):
    X = np.genfromtxt(ifile, skip_header=5)
    xcoor = X[:, 0]
    ycoor = X[:, 1]
    zcoor = X[:, 2]
    return xcoor, ycoor, zcoor


def write_files(proc, mask_noise, normal_x_noise, normal_y_noise,
                normal_z_noise, outdir):
    _write(mask_noise,
           os.path.join(outdir, 'proc{:06}_mask_noise.bin'.format(proc)))
    _write(normal_x_noise,
           os.path.join(outdir, 'proc{:06}_normal_x_noise.bin'.format(proc)))
    _write(normal_y_noise,
           os.path.join(outdir, 'proc{:06}_normal_y_noise.bin'.format(proc)))
    _write(normal_z_noise,
           os.path.join(outdir, 'proc{:06}_normal_z_noise.bin'.format(proc)))
    return


def _write(x, filename):
    n = np.array([4 * len(x)], dtype='int32')
    x = np.array(x, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        x.tofile(file)
        n.tofile(file)


if __name__ == '__main__':
    fname = sys.argv[1]

    try:
        with open(fname, 'r') as _file:
            par = yaml.load(_file)
    except IOError:
        print('IOError: parfile not found.')

    if par['PLOT_MASK']:
        all_gll_x = np.array([])
        all_gll_y = np.array([])
        all_mask_noise = np.array([])

    if par['ADD_RANDOM_DISTRIBUTION']:
        dx = max((par['XMAX']-par['XMIN']),
                 (par['YMAX']-par['YMIN'])) / par['NX']

        random_field = gaussian_random_field(par['NX'], dx, par['WAVELENGHT'])
        random_field[np.where(random_field < 0.0)] = 0.0
        random_field /= np.max(random_field)
        random_field *= par['AMP_FACTOR']

        x2 = np.linspace(par['XMIN'], par['XMAX'], par['NX'], endpoint=True)
        y2 = np.linspace(par['YMIN'], par['YMAX'], par['NX'], endpoint=True)

        interpolator = RegularGridInterpolator((x2, y2), random_field,
                                               method='linear')

    for p in range(0, par['NPROC']):
        if par['FREE_SURFACE_PROC'] and p not in par['FREE_SURFACE_PROC']:
            if par['WRITE_FILES']:
                print('Writting dummy files for process {}'.format(p))
                write_files(p, [0.0], [0.0], [0.0], [0.0],
                            par['DATABASES_MPI'])
            continue

        print('Generating files for process {}'.format(p))
        # read gll coordinates
        free_surface_file = os.path.join(par['DATABASES_MPI'],
                                         'proc{:06}_free_surface.vtk'.format(
                                             p))

        gll_x, gll_y, gll_z = read_gll_coordinates(free_surface_file)

        # initialize arrays
        proc_ngll_surface = gll_x.size
        mask_noise = np.zeros(proc_ngll_surface, dtype='float32')
        normal_x_noise = np.zeros(proc_ngll_surface, dtype='float32')
        normal_y_noise = np.zeros(proc_ngll_surface, dtype='float32')
        normal_z_noise = np.zeros(proc_ngll_surface, dtype='float32')

        # set noise direction
        normal_x_noise[:] = par['XDIR']
        normal_y_noise[:] = par['YDIR']
        normal_z_noise[:] = par['ZDIR']

        # set noise distribution
        for d in par['DISTRIBUTIONS']:
            if d['type'] == 'uniform':
                mask_noise = uniform_distribution(d, mask_noise)
            elif d['type'] == 'ocean':
                mask_noise = ocean_distribution(d, mask_noise, gll_z)
            elif d['type'] == 'gaussian':
                mask_noise = gaussian_distribution(d, mask_noise, gll_x, gll_y)
            else:
                print('Undefined noise distribution.')

        if par['ADD_RANDOM_DISTRIBUTION']:
            points = np.hstack((gll_x.reshape(-1, 1), gll_y.reshape(-1, 1)))
            mask_noise += interpolator(points)

            if par['DONT_PERTURB']:
                rcoor = np.sqrt(np.square(gll_x - par['X0']) +
                                np.square(gll_y - par['Y0']))

                mute = (rcoor - par['R1']) / (par['R2'] - par['R1'])

                idx = np.argwhere(rcoor <= par['R1'])
                mute[idx] = 0.0

                idx = np.argwhere(rcoor >= par['R2'])
                mute[idx] = 1.0

                mask_noise *= mute

        if par['WRITE_FILES']:
            write_files(p, mask_noise, normal_x_noise, normal_y_noise,
                        normal_z_noise, par['OUTPUT_PATH'])

        if par['PLOT_MASK']:
            all_gll_x = np.append(all_gll_x, gll_x)
            all_gll_y = np.append(all_gll_y, gll_y)
            all_mask_noise = np.append(all_mask_noise, mask_noise)

    if par['PLOT_MASK']:
        fig, ax = plt.subplots()
        im = ax.scatter(all_gll_x, all_gll_y, c=all_mask_noise)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_title('Noise distribution mask')
        cbar = plt.colorbar(im)
        cbar.set_label('Weight')
        plt.show()
        plt.close()
