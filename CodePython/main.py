#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:37:59 2020

@author: quenot
"""

import sys
import datetime
import time
sys.path.append('CodePython/InputOutput/')
import os
from InputOutput.pagailleIO import saveEdf
from Experiment import Experiment
# from ExperimentMultiprocess import Experiment
import numpy as np


if __name__ == "__main__":
    """main of the simulation code.

    Notes:
        Set the parameters below and parameters in .xml files then launch
    """
    time0=time.time() #timer for computation
    exp_dict={}
    
    ## PARAMETERS TO SET
    # Define experiment 
    exp_dict['experimentName']="SIMAP_SpheresInTube"
    # Output filepath to store the result images
    exp_dict['filepath']='../Results/SIMAP_SpheresInTube/'
    # Define algorithm parameters
    exp_dict['sampleSampling']=2 # MUST BE AN INTEGER
    exp_dict['margin']=10 #with Fresnel there might be an aliasing issue so we need to extend very slightly the image for calculations
    save=True
    exp_dict['simulation_type']="RayT" #"Fresnel" or "RayT" 

   
    #************************************************************************
    now=datetime.datetime.now()
    exp_dict['expID']=now.strftime("%Y%m%d-%H%M%S") #define experiment ID
    
    PropagImage=[]
    Geometry=[]
    
    experiment=Experiment(exp_dict) 
    
    print("\n\nINITIALIZING EXPERIMENT PARAMETERS AND GEOMETRIES")
    print("\n\n*************************")
    if exp_dict['simulation_type']=="Fresnel":
        PropagImageTmp, White=experiment.computeSampleAndReferenceImages()
    elif exp_dict['simulation_type']=="RayT":
        PropagImageTmp, White=experiment.computeSampleAndReferenceImagesRT()
    else:
        raise Exception("simulation Type not defined: ", exp_dict['simulation_type'])
    Nbin=len(PropagImageTmp)
    
    expPathEn=[]
    if exp_dict['simulation_type']=="Fresnel":
        expImagesFilePath=exp_dict['filepath']+'Fresnel_'+str(exp_dict['expID'])+'/'
    if exp_dict['simulation_type']=="RayT":
        expImagesFilePath=exp_dict['filepath']+'RayTracing_'+str(exp_dict['expID'])+'/'
    os.mkdir(expImagesFilePath)
    thresholds=experiment.myDetector.myBinsThersholds.copy()
    thresholds.insert(0,experiment.mySource.mySpectrum[0][0])
    for ibin in range(Nbin):
        binstart='%2.2d'%thresholds[ibin]
        binend='%2.2d'%thresholds[ibin+1]
        expPathEn.append(f'{expImagesFilePath}{binstart}_{binend}kev/')
        if len(thresholds)-1==1:
            expPathEn=[expImagesFilePath]
        else:
            os.mkdir(expPathEn[ibin])
                
    for ibin in range(Nbin):
        saveEdf(PropagImageTmp[ibin], expPathEn[ibin]+'PropagImage_'+str(exp_dict['expID'])+'.edf')
        saveEdf(White[ibin], expPathEn[ibin]+'White_'+str(exp_dict['expID'])+'.edf')

    experiment.saveAllParameters(time0,exp_dict)
    
    print("\nfini")
    



        