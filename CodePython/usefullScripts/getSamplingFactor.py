#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 11:21:43 2021
This code calculates the sampling requirements for object simulation

Inspired from Häggmark, I., Shaker, K., & Hertz, H. M. (2020). In Silico Phase-Contrast X-Ray Imaging of Anthropomorphic Voxel-Based Phantoms. IEEE Transactions on Medical Imaging, 40(2), 539-548.
@author: quenot
"""
import numpy as np

def kevToLambda(energyInKev):
    energy = energyInKev * 1e3
    waveLengthInNanometer = 1240. / energy
    return waveLengthInNanometer * 1e-9

def is_overSampling_ok(exp_dict, pixel_size, energy):
    if exp_dict['simulation_type']=="Fresnel":
        magnification=(exp_dict['distSourceToMembrane']+exp_dict['distMembraneToObject']+exp_dict['distObjectToDetector'])/(exp_dict['distSourceToMembrane']+exp_dict['distMembraneToObject'])
        lambda_val=kevToLambda(energy)
        min_dx=np.sqrt(lambda_val*exp_dict['distObjectToDetector']/magnification)/2
        pix_size_in_sample_plane=pixel_size/magnification
        min_oversampling=np.ceil(pix_size_in_sample_plane/min_dx*1e-6)
        if min_oversampling>exp_dict['overSampling']:
            print(f'/!\/!\ OVERSAMPLING FACTOR < MIN OVERSAMPLING FOR FRESNEL MODEL: {exp_dict["overSampling"]} < {min_oversampling}')
    return min_oversampling
    

if __name__ == "__main__":
    # PARAMETERS
    exp_dict={}
    exp_dict['Energy'] = 22 #keV
    exp_dict['overSampling']=2 #The function need a value, if you do not know it yet put a random integer or leave 2
    exp_dict['distSourceToMembrane']=0.5 #in m
    exp_dict['distMembraneToObject']=0 #in m
    exp_dict['distObjectToDetector']=1 #in m
    exp_dict['simulation_type']="Fresnel"
    detector_pixel_size=50 #in um
    
    
    min_oversamp=is_overSampling_ok(exp_dict, detector_pixel_size, exp_dict['Energy'])
    
    print('Min oversampling factor:', min_oversamp)
