import matplotlib.pyplot as plt
import numpy as np


def fft_idx_gen(size):
    '''
    generate fft indexes
    '''

    a = np.arange(0, size//2+1)
    b = np.arange(1, size//2)
    b = b[::-1] * -1
    kidx = np.append(a, b)
    return kidx


def gaussian_spectrum(kx, ky, corr_len):
    mean = 0.0
    std = 1.0 / (2.0 * corr_len)

    k = np.sqrt(kx**2 + ky**2)
    amp = np.exp(-(k-mean)**2/(2*std**2))

    return amp


def gaussian_random_field(nx, dx, corr_len):
    kres = 1.0 / (nx * dx)  # wavenumber resolution
    kaxis = fft_idx_gen(nx) * kres  # wavenumber axis

    noise = np.random.normal(size=(nx, nx))
    noise_fft = np.fft.fftn(noise)

    amplitude = np.zeros((nx, nx))

    for i, kx in enumerate(kaxis):
        for j, ky in enumerate(kaxis):
            amplitude[i, j] = gaussian_spectrum(kx, ky, corr_len)

    plt.plot(kaxis, amplitude[i, :], 'o')
    plt.title('Wavenumber spectrum')
    plt.show()

    random_field = np.fft.ifftn(noise_fft * amplitude)
    return np.real(random_field)


if __name__ == '__main__':
    corr_len = 50.0
    nx = 100
    dx = 1.0

    perturbations = gaussian_random_field(nx, dx, corr_len)

    plt.imshow(perturbations, extent=[0, nx*dx, 0, nx*dx])
    plt.colorbar()
    plt.show()
