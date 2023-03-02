import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def fft_idx_gen(size):
    '''
    generate fft indexes
    '''

    a = np.arange(0, size//2+1)
    b = np.arange(1, size//2)
    b = b[::-1] * -1
    kaxis = np.append(a, b)
    return kaxis


def correlator_spectrum(kx, ky, kz, alpha):
    '''
    power spectrum of the form 1/k^alpha
    '''
    if kx == 0 and ky == 0 and kz ==0:
        amp = 0.0
    else:
        k = np.sqrt(kx**2 + ky**2 + kz**2)
        power = k**alpha
        amp = np.sqrt(power)
    return amp


def gaussian_random_field(alpha, size=100):
    kidx = fft_idx_gen(size)

    noise = np.random.normal(loc=0.0,
                             scale=1.0,
                             size=(size, size, size))

    noise_fft = np.fft.fftn(noise)

    amplitude = np.zeros((size, size, size))

    for i, kx in enumerate(kidx):
        for j, ky in enumerate(kidx):
            for k, kz in enumerate(kidx):
                amplitude[i, j, k] = correlator_spectrum(kx, ky, kz, alpha)

    random_field = np.real(np.fft.ifftn(noise_fft * amplitude))
    return random_field


if __name__ == '__main__':
    # USER PARAMETERS
    # ===============
    tomo_file = './tomography_model.xyz'

    # random field 
    pfactor = 0.1
    alpha = -3
    nel = 128

    plot_figures = True

    # DO NOT EDIT BELOW THIS LINE
    # ===========================

    # generate perturbations
    perturbations = gaussian_random_field(alpha, size=nel)

    # scale perturbations between -pfactor and pfactor
    perturbations = perturbations / np.max(np.abs(perturbations))
    perturbations *= pfactor

    # read model to perturb
    header = np.loadtxt(tomo_file, skiprows=2, max_rows=1, dtype='int')
    model = np.loadtxt(tomo_file, skiprows=4, dtype='float32')

    nx = header[0]
    ny = header[1]
    nz = header[2]

    x = model[:, 0].reshape(-1, 1)
    y = model[:, 1].reshape(-1, 1)
    z = model[:, 2].reshape(-1, 1)

    vp = model[:, 3].copy()
    vs = model[:, 4].copy()
    rho = model[:, 5].copy()

    # generate interpolator object
    x2 = np.linspace(x.min(), x.max(), nel, endpoint=True)
    y2 = np.linspace(y.min(), y.max(), nel, endpoint=True)
    z2 = np.linspace(z.min(), z.max(), nel, endpoint=True)

    interpolator = RegularGridInterpolator((x2, y2, z2), perturbations,
                                           method='nearest')

    # perturb model
    points = np.hstack((x, y, z))
    dm = interpolator(points)

    model[:, 3] = vp * (1.0 + dm)
    model[:, 4] = vs * (1.0 + dm)
    model[:, 5] = rho * (1.0 + dm)

    print('Min perturbation: {} %'.format(dm.min()))
    print('Max perturbation: {} %\n'.format(dm.max()))

    print('Initial Vp min: {}'.format(vp.min()))
    print('Initial Vp max: {}\n'.format(vp.max()))

    print('Perturbed Vp min: {}'.format(model[:, 3].min()))
    print('Perturbed Vp max: {}\n'.format(model[:, 3].max()))

    # save perturbed model
    np.savetxt('tomography_model.xyz.perturbed', model, fmt='%1.3f')

    # figures
    if plot_figures:
        yidx = 0

        plt.imshow(perturbations[:, yidx, :].T, origin='lower')
        plt.xlabel('X')
        plt.ylabel('Z')
        plt.title('Perturbations: vertical slice')
        plt.colorbar()
        plt.show()

        vp = np.reshape(vp, (nx, ny, nz), order='F')

        plt.imshow(vp[:, yidx, :].T, origin='lower')
        plt.xlabel('X')
        plt.ylabel('Z')
        plt.title('Initial Vp: vertical slice')
        plt.colorbar()
        plt.show()

        vp2 = np.reshape(model[:, 3], (nx, ny, nz), order='F')

        plt.imshow(vp2[:, yidx, :].T, origin='lower')
        plt.xlabel('X')
        plt.ylabel('Z')
        plt.title('Perturbed Vp: vertical slice')
        plt.colorbar()
        plt.show()
