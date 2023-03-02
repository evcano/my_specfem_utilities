import numpy as np
import os
from mpi4py import MPI
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt


def fft_idx_gen(size):
    a = np.arange(0, size//2+1)
    b = np.arange(1, size//2)
    b = b[::-1] * -1
    kaxis = np.append(a, b)
    return kaxis


def gaussian_spectrum(kx, ky, kz, std):
    mean = 0.0
    k = np.sqrt(kx**2 + ky**2 + kz**2)
    amp = np.exp(-(k-mean)**2 / (2*std**2))
    return amp


def gaussian_random_field(nx, dx, wavelength):
    std = 1.0 / (2.0 * wavelength)

    kres = 1.0 / (nx * dx)
    kidx = fft_idx_gen(nx)
    kaxis = kidx * kres

    noise = np.random.normal(loc=0.0, scale=1.0, size=(nx, nx, nx))
    noise_fft = np.fft.fftn(noise)

    amplitude = np.zeros((nx, nx, nx))

    for h, kx in enumerate(kaxis):
        for i, ky in enumerate(kaxis):
            for j, kz in enumerate(kaxis):
                amplitude[h, i, j] = gaussian_spectrum(kx, ky, kz, std)

    random_field = np.real(np.fft.ifftn(noise_fft * amplitude))
    return random_field


def _write(x, filename):
    n = np.array([4 * len(x)], dtype='int32')
    x = np.array(x, dtype='float32')
    with open(filename, 'wb') as file:
        n.tofile(file)
        x.tofile(file)
        n.tofile(file)


def _read(filename, dtype='float32'):
    nbytes = os.path.getsize(filename)
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype='int32', count=1)[0]

        if n == nbytes-8:
            file.seek(4)
            data = np.fromfile(file, dtype=dtype)
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype=dtype)
            return data


if __name__ == '__main__':
    # MPI
    comm = MPI.COMM_WORLD
    proc = comm.Get_rank()
    nproc = comm.Get_size()

    # USER PARAMETERS
    # ===============
    databases_mpi = '../../DATABASES_MPI'
    outdir = '.'
    model_parameters = ['vp', 'vs']

    # mesh info
    xmin = 283058.188
    xmax = 1003078.19
    ymin = 1800111.12
    ymax = 2520131.00
    zmin = -240000.0
    zmax = 0.0

    # random gaussian field
    wavelength = 22000.0
    pstd = 0.20  # perturbations std
    pmax = 0.25  # maximum perturbation
    nx = 128

    # dont perturb area
    dont_perturb = True
    R1 = 80000.0
    R2 = 120000.0
    R3 = 200000.0
    a_x0 = 643068.1804806339
    a_y0 = 2160121.0842638565
    a_zmin = zmin
    a_zmax = zmax

    # DO NOT EDIT BELOW THIS LINE
    # ===========================

    if proc == 0:
        dx = max((xmax-xmin), (ymax-ymin), (zmax-zmin)) / nx

        print('generating gaussian perturbations ...')
        perturbations = gaussian_random_field(nx, dx, wavelength)

        # scale perturbations between -pfactor and pfactor
        perturbations = perturbations / np.max(np.abs(perturbations))
        perturbations *= pstd * 3.0
        perturbations[np.where(perturbations > pmax)] = pmax
        perturbations[np.where(perturbations < -pmax)] = -pmax

        plt.hist(perturbations.reshape(-1,1))
        plt.show()

        plt.imshow(perturbations[:, :, 0])
        plt.show()

        # generate interpolator instance
        x2 = np.linspace(xmin, xmax, nx, endpoint=True)
        y2 = np.linspace(ymin, ymax, nx, endpoint=True)
        z2 = np.linspace(zmin, zmax, nx, endpoint=True)

        print('generating perturbations interpolator ...')
        interpolator = RegularGridInterpolator((x2, y2, z2), perturbations,
                                               method='linear')
    else:
        interpolator = None

    interpolator = comm.bcast(interpolator, root=0)

    print(' perturbing model ...')

    # read gll coordinates
    ibool_file = 'proc{:06}_ibool.bin'.format(proc)
    xcoor_file = 'proc{:06}_x.bin'.format(proc)
    ycoor_file = 'proc{:06}_y.bin'.format(proc)
    zcoor_file = 'proc{:06}_z.bin'.format(proc)

    ibool = _read(os.path.join(databases_mpi, ibool_file), dtype='int32')
    xcoor = _read(os.path.join(databases_mpi, xcoor_file), dtype='float32')
    ycoor = _read(os.path.join(databases_mpi, ycoor_file), dtype='float32')
    zcoor = _read(os.path.join(databases_mpi, zcoor_file), dtype='float32')

    # interpolate perturbations on the gll points
    gll_points = np.hstack((xcoor.reshape(-1, 1),
                            ycoor.reshape(-1, 1),
                            zcoor.reshape(-1, 1)))

    dm_all = interpolator(gll_points)

    # the model in these gll points wil not be perturbed
    if dont_perturb:
        rcoor = np.sqrt(np.square(xcoor - a_x0) + np.square(ycoor - a_y0))
        
        idx_local = np.argwhere(rcoor <= R1)
        idx_transition1 = np.argwhere((rcoor <= R2) & (rcoor > R1))
        idx_transition2 = np.argwhere((rcoor < R3) & (rcoor > R2))
        idx_background = np.argwhere(rcoor > R3)
        
        mute = np.zeros(rcoor.shape)
        mute[idx_local] = 1.0 
        tmp = abs(R2 - rcoor[idx_transition1])
        mute[idx_transition1] = (tmp - np.min(tmp)) / (np.max(tmp) - np.min(tmp))
        tmp = rcoor[idx_transition2]
        mute[idx_transition2] = (tmp - np.min(tmp)) / (np.max(tmp) - np.min(tmp))
        mute[idx_background] = 1.0
        plt.plot(rcoor, mute, 'o')
        plt.show()
        
        dm_all *= mute
        
        dm_background = dm_all.copy()
        idx_local = np.argwhere(rcoor <= R2)
        dm_background[idx_local] = 0.0
        
        dm_local = dm_all.copy()
        idx_background = np.argwhere(rcoor > R2)
        dm_local[idx_background] = 0.0

    # read and perturb  model
    ibool -= 1

    # perturb all model
    for mpar in model_parameters:
        m_file = 'proc{:06}_{}.bin'.format(proc, mpar)
        m = _read(os.path.join(databases_mpi, m_file), dtype='float32')

        m_all = m * (1.0 + dm_all[ibool])

        outdir2 = os.path.join(outdir, 'all')
        if not os.path.isdir(outdir2):
            os.makedirs(outdir2)

        _write(m_all, os.path.join(outdir2, f'{m_file}'))

        if dont_perturb:
            m_background = m * (1.0 + dm_background[ibool])

            outdir2 = os.path.join(outdir, 'background')
            if not os.path.isdir(outdir2):
                os.makedirs(outdir2)

            _write(m_background, os.path.join(outdir2, f'{m_file}'))

            m_local = m * (1.0 + dm_local[ibool])

            outdir2 = os.path.join(outdir, 'local')
            if not os.path.isdir(outdir2):
                os.makedirs(outdir2)

            _write(m_local, os.path.join(outdir2, f'{m_file}'))

    print('proc {} done'.format(proc))
