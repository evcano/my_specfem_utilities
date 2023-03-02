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


def read_gll_coordinates(ifile):
    X = np.genfromtxt(ifile, skip_header=5)
    xcoor = X[:, 0]
    ycoor = X[:, 1]
    zcoor = X[:, 2]
    return xcoor, ycoor, zcoor


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


DATABASES_MPI = './DATABASES_MPI'
NSPEC_SURFACE_EXT_MESH = 1024
NGLLSQUARE = 25
NPROC = 8

if __name__ == '__main__':
    ngll_surface = NSPEC_SURFACE_EXT_MESH * NGLLSQUARE
    proc_ngll_surface = ngll_surface // NPROC

    all_gll_x = np.array([])
    all_gll_y = np.array([])
    all_mask_noise = np.array([])

    for p in range(0, NPROC):

        # read gll coordinates
        free_surface_file = os.path.join(DATABASES_MPI,
                                         'proc{:06}_free_surface.vtk'.format(p))

        gll_x, gll_y, gll_z = read_gll_coordinates(free_surface_file)

        all_gll_x = np.append(all_gll_x, gll_x)
        all_gll_y = np.append(all_gll_y, gll_y)

        # read noise mask
        mask_noise = np.zeros(proc_ngll_surface, dtype='float32')

        mask_noise_file = os.path.join(DATABASES_MPI,
                                       'proc{:06}_mask_noise.bin'.format(p))
        mask_noise = _read(mask_noise_file)

        all_mask_noise = np.append(all_mask_noise, mask_noise)

    # figure
    fig, ax = plt.subplots()
    im = ax.scatter(all_gll_x, all_gll_y, c=all_mask_noise)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title('Noise distribution mask')
    cbar = plt.colorbar(im)
    cbar.set_label('Weight')
    plt.show()
    plt.close()
