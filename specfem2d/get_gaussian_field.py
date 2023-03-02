#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import pickle
import sys
import yaml
from gstools import SRF, Gaussian
from scipy.interpolate import RegularGridInterpolator


if __name__ == '__main__':
    parfile = sys.argv[1]

    try:
        with open(parfile, 'r') as _file:
            par = yaml.safe_load(_file)
    except IOError:
        print('IOError: parfile not found.')

    # DO NOT EDIT BELOW THIS LINE
    model = Gaussian(dim=2, var=par['DSTD']**2,
                     len_scale=par['LEN_SCALE'])

    srf = SRF(model, mean=par['DMEAN'])

    x = np.linspace(par['XMIN'], par['XMAX'], par['NX'], endpoint=True)
    z = np.linspace(par['ZMIN'], par['ZMAX'], par['NZ'], endpoint=True)

    perturbations = srf.structured([x, z])
    print(np.mean(perturbations))
    print(perturbations.min(), perturbations.max())

    print('generating perturbations interpolator ...')
    interpolator = RegularGridInterpolator((x, z), perturbations, method='linear')

    # save iterpolator object
    output_file = os.path.join(par['OUTPUT_DIR'], par['NAME'])

    with open(f'{output_file}.pkl', 'wb') as _file:
        pickle.dump(interpolator, _file)

    # plot field
    cmap = cm.get_cmap('seismic')
    plt.imshow(perturbations, cmap=cmap, vmin=-1.0, vmax=1.0)
    plt.colorbar()
    plt.show()
