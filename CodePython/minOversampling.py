import numpy as np

def kevToLambda(energyInKev):
    energy = energyInKev * 1e3
    waveLengthInNanometer = 1240. / energy
    return waveLengthInNanometer * 1e-9

def minOversampling(sample_detector, source_sample, energy, pixel_size):
    magnification = (sample_detector+source_sample)/source_sample
    wavelength = kevToLambda(energy)
    min_dx = np.sqrt(wavelength*sample_detector/magnification)/2
    min_oversampling = int(-((pixel_size/min_dx/magnification*1e-6)//-1))
    return min_oversampling
