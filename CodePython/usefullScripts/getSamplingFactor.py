#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 11:21:43 2021
This code calculates the sampling requirements for object simulation

Inspired from Häggmark, I., Shaker, K., & Hertz, H. M. (2020). In Silico Phase-Contrast X-Ray Imaging of Anthropomorphic Voxel-Based Phantoms. IEEE Transactions on Medical Imaging, 40(2), 539-548.
@author: quenot
"""
import numpy as np
# PARAMETERS
energy = 39.5 #keV
dist_source_sample=0.5 #in m
dist_sample_detector=1 #in m
detector_pixel_size=50 #in um


def kevToLambda(energyInKev):
    energy = energyInKev * 1e3
    waveLengthInNanometer = 1240. / energy
    return waveLengthInNanometer * 1e-9

magnification=(dist_sample_detector+dist_source_sample)/dist_source_sample
lambda_val=kevToLambda(energy)

min_dx=np.sqrt(lambda_val*dist_sample_detector/magnification)/2

print("minimal sample sampling size", min_dx)

pix_size_in_sample_plane=detector_pixel_size/magnification
min_oversampling=np.ceil(pix_size_in_sample_plane/min_dx*1e-6)

print("minimal oversampling", min_oversampling)

