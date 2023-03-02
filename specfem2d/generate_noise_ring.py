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

# FOR SPECFEM2D !!!

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import yaml
from scipy.interpolate import griddata, interp1d


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


def _read(filename):
    """
    Reads Fortran style binary data into numpy array
    """
    nbytes = os.path.getsize(filename)
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype='int32', count=1)[0]

        if n == nbytes-8:
            file.seek(4)
            data = np.fromfile(file, dtype='float32')
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype='float32')
            return data


if __name__ == '__main__':
    fname = sys.argv[1]

    try:
        with open(fname, 'r') as _file:
            par = yaml.load(_file)
    except IOError:
        print('IOError: parfile not found.')

    ngll_boundary = par['NSPEC_THETA'] * 5
    ngll_r = par['NSPEC_R'] * 5

    proc_ngll = []
    all_gll_x = np.array([])
    all_gll_z = np.array([])
    all_mask_noise = np.array([])

    # read gll points coordinates
    for p in range(0, par['NPROC']):
        mesh_glob = np.loadtxt(os.path.join(par['DATADIR'],
                                            'mesh_glob{:06}'.format(p)))

        gll_x = mesh_glob[:, 0]
        gll_z = mesh_glob[:, 1]

        proc_ngll.append(gll_x.size)

        all_gll_x = np.append(all_gll_x, gll_x)
        all_gll_z = np.append(all_gll_z, gll_z)

    ngll = all_gll_x.size

    # compute radial coordinate wrt the mesh geometrical center
    all_gll_r = np.sqrt((all_gll_x - par['MESH_X0'])**2 +
                        (all_gll_z - par['MESH_Z0'])**2)

    rmax = np.max(all_gll_r)

    # determine the gll points of the noise-source ring
    ring_limit = rmax - par['RING_THICKNESS']
    ring_idx = np.argwhere(all_gll_r >= ring_limit)

    # compute radial and angular coordinates wrt the integral reference point
    all_gll_r = np.sqrt((all_gll_x - par['REF_X0'])**2 +
                        (all_gll_z - par['REF_Z0'])**2)

    all_gll_theta = (np.arctan2(all_gll_z - par['REF_Z0'], all_gll_x -
                                par['REF_X0']) + 2*np.pi) % (2*np.pi)

    rmax = np.max(all_gll_r)

    # obtain noise ring by integrating an existing noise mask
    if par['INTEGRATE']:
        # read noise mask
        for p in range(0, par['NPROC']):
            mask_noise = _read(os.path.join(par['NOISE_DISTRIBUTION'],
                                            'proc{:06}_mask_noise.bin'.format(
                                                p)))

            all_mask_noise = np.append(all_mask_noise, mask_noise)

        # generate a polar grid centered at the integral reference point
        dr = rmax / ngll_r
        dtheta = (2*np.pi) / ngll_boundary

        ax_r = np.arange(0.0, rmax + dr, dr)
        ax_theta = np.arange(0.0, 2*np.pi + dtheta, dtheta)

        R, T = np.meshgrid(ax_r, ax_theta)
        X = R * np.cos(T) + par['REF_X0']
        Z = R * np.sin(T) + par['REF_Z0']

        # interpolate noise mask into the polar grid
        points = np.zeros((ngll, 2))
        points[:, 0] = np.squeeze(all_gll_x)
        points[:, 1] = np.squeeze(all_gll_z)

        mask_noise_intp = griddata(points, all_mask_noise, (X, Z),
                                   method='linear', fill_value=np.nan)

        # azimuthal integration to obtain the noise-source ring
        r = np.zeros(ax_r.size)
        r[1:] = 0.5 * (ax_r[0:-1] + ax_r[1:])
        da = r * dr * dtheta

        noise_ring = np.zeros(ax_theta.size)

        for i, theta in enumerate(ax_theta):
            idx = np.where(T == theta)
            noise_ring[i] = np.nansum(mask_noise_intp[idx] * da)

    # read noise ring from file
    else:
        tmp = np.loadtxt(par['NOISE_DISTRIBUTION'])
        ax_theta = np.deg2rad(tmp[:, 0])
        noise_ring = tmp[:, 1]

    np.savetxt('noise_ring.output', np.array([np.rad2deg(ax_theta),
                                              noise_ring]).T)

    # interpolate noise ring to gll points
    f = interp1d(ax_theta, noise_ring)
    all_mask_noise_boundary = f(all_gll_theta[ring_idx])

    # get noise mask with values only at the ring
    all_mask_noise_ring = np.zeros(ngll)
    all_mask_noise_ring[ring_idx] = all_mask_noise_boundary

    # save files
    if par['WRITE_FILES']:
        idx1 = 0
        idx2 = 0

        for p in range(0, par['NPROC']):
            idx1 = idx2
            idx2 += proc_ngll[p]
            mask_noise = all_mask_noise_ring[idx1:idx2]
            _write(mask_noise,
                   os.path.join(par['OUTPUT_PATH'], 'proc{:06}_mask_noise.bin'.format(p)))

    # figures
    if par['PLOT_MASK']:
        if par['INTEGRATE']:
            fig, ax = plt.subplots()
            im = ax.tripcolor(all_gll_x, all_gll_z, all_mask_noise)
            ax.set_xlabel('x [m]')
            ax.set_ylabel('z [m]')
            ax.set_title('Input noise distribution mask')
            cbar = plt.colorbar(im)
            cbar.set_label('Weight')
            plt.show()
            plt.close()

        fig, ax = plt.subplots()
        ax.scatter(np.rad2deg(ax_theta), noise_ring)
        ax.set_xlabel('Angle in degrees')
        ax.set_ylabel('Weight')
        ax.set_title('Noise directionality')
        plt.show()
        plt.close()

        fig, ax = plt.subplots()
        im = ax.tripcolor(all_gll_x, all_gll_z, all_mask_noise_ring)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('z [m]')
        ax.set_title('Output noise distribution mask')
        cbar = plt.colorbar(im)
        cbar.set_label('Weight')
        plt.show()
        plt.close()
