"""Calculate the minimum sampling requirement for Fresnel simulation
Inspired from Häggmark, I., Shaker, K., & Hertz, H. M. (2020).
In Silico Phase-Contrast X-Ray Imaging of Anthropomorphic Voxel-Based Phantoms.
IEEE Transactions on Medical Imaging, 40(2), 539-548.
"""
import numpy as np


def kev_lambda(energy_kev):
    """Convert between energy in keV and wavelength in nm and vice-versa

    Parameters
    ----------
    energy_kev : float
        Energy in keV or wavelength in nm

    Returns
    -------
    float
        Wavelength in nm or energy in keV
    """
    energy = energy_kev * 1e3
    wavelength_nm = 1240. / energy
    return wavelength_nm * 1e-9


def min_oversampling(sample_detector, source_sample, energy, pixel_size):
    """Calculate minimum oversampling requirement for Fresnel simulations

    Parameters
    ----------
    sample_detector : float
        Sample to detector distance in m
    source_sample : float
        Source to sample distance in m
    energy : float
        Energy in keV
    pixel_size : float
        Pixel size in μm

    Returns
    -------
    int
        The minimum oversampling factor for the given Fresnel simulation
    """
    magnification = (sample_detector+source_sample)/source_sample
    wavelength = kev_lambda(energy)
    min_dx = np.sqrt(wavelength*sample_detector/magnification)/2
    # Ceiling operator
    min_sampling = int(-((pixel_size/min_dx/magnification*1e-6)//-1))
    return min_sampling
