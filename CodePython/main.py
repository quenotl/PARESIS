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
from InputOutput.pagailleIO import save_image
from Experiment import Experiment
# from ExperimentMultiprocess import Experiment
import numpy as np


if __name__ == "__main__":
    """main of the simulation code.

    Notes:
        Set the parameters below and parameters in .xml files then run
    """
    time0=time.time() #timer for computation
    exp_dict={}
    
    ## PARAMETERS TO SET
    # Define experiment 
    exp_dict['experimentName']="Fil_Nylon_ID17"
    # Output filepath to store the result images
    exp_dict['filepath']='../Results/Fil_Nylon_ID17/'
    # Define algorithm parameters
    exp_dict['overSampling']=2 # MUST BE AN INTEGER - >2 for ray-tracing model and even more for Fresnel (cf usefullScripts/getSamplingFactor.py)
    exp_dict['nbExpPoints']=1 #number of pair of acquisitions (Ir, Is) simulated with different positions of the membrane
    save=True
    saving_format='.tif' #.tif or .edf
    exp_dict['simulation_type']="RayT" #"Fresnel" or "RayT" 

   
    #************************************************************************
    #**********START OF CALCULATIONS*****************************************
    #************************************************************************
    
    now=datetime.datetime.now()
    exp_dict['expID']=now.strftime("%Y%m%d-%H%M%S") #define experiment ID
    
    SampleImage=[]
    ReferenceImage=[]
    PropagImage=[]
    AbsImage=[]
    SubImage=[]
    Geometry=[]
    
    print("\n\nINITIALIZING EXPERIMENT PARAMETERS AND GEOMETRIES")
    print("*************************")
    experiment=Experiment(exp_dict) 
    
    
    print("\nImages calculation")
    print("*************************")
    for pointNum in range(exp_dict['nbExpPoints']):
        experiment.myMembrane.myGeometry=[]
        experiment.myMembrane.getMyGeometry(experiment.exp_dict['studyDimensions'],experiment.myMembrane.membranePixelSize,experiment.exp_dict['overSampling'], pointNum, exp_dict['nbExpPoints'])

        print("\nCalculations point",pointNum)
        if exp_dict['simulation_type']=="Fresnel":
            SampleImageTmp, ReferenceImageTmp,PropagImageTmp, White=experiment.computeSampleAndReferenceImages_Fresnel(pointNum)
        elif exp_dict['simulation_type']=="RayT":
            SampleImageTmp, ReferenceImageTmp,PropagImageTmp, White, Dx, Dy, Df=experiment.computeSampleAndReferenceImages_RT(pointNum)
        else:
            raise Exception("simulation Type not defined: ", exp_dict['simulation_type'])
        Nbin=len(SampleImageTmp)
        
        if pointNum==0:
            expPathEn=[]
            if exp_dict['simulation_type']=="Fresnel":
                expImagesFilePath=exp_dict['filepath']+'Fresnel_'+str(exp_dict['expID'])+'/'
            if exp_dict['simulation_type']=="RayT":
                expImagesFilePath=exp_dict['filepath']+'RayTracing_'+str(exp_dict['expID'])+'/'
            os.mkdir(expImagesFilePath)
            os.mkdir(expImagesFilePath+'membraneThickness/')
            thresholds=experiment.myDetector.det_param['myBinsThersholds'].copy()
            thresholds.insert(0,experiment.mySource.mySpectrum[0][0])
            for ibin in range(Nbin):
                binstart='%2.2d'%thresholds[ibin]
                binend='%2.2d'%thresholds[ibin+1]
                expPathEn.append(f'{expImagesFilePath}{binstart}_{binend}kev/')
                if len(thresholds)-1==1:
                    expPathEn=[expImagesFilePath]
                else:
                    os.mkdir(expPathEn[ibin])
                os.mkdir(expPathEn[ibin]+'ref/')
                os.mkdir(expPathEn[ibin]+'sample/')
                os.mkdir(expPathEn[ibin]+'propag/')
            
        txtPoint = '%2.2d' % pointNum
        save_image(experiment.myMembrane.myGeometry[0], expImagesFilePath+'membraneThickness/'+exp_dict['experimentName']+'_sampling'+str(exp_dict['overSampling'])+'_'+str(pointNum)+saving_format)
        
        if exp_dict['simulation_type']=="RayT":
            save_image(Df, expImagesFilePath+"DF"+saving_format)
            
        for ibin in range(Nbin):
            save_image(SampleImageTmp[ibin], expPathEn[ibin]+'sample/sampleImage_'+str(exp_dict['expID'])+'_'+txtPoint+saving_format)
            save_image(ReferenceImageTmp[ibin], expPathEn[ibin]+'ref/ReferenceImage_'+str(exp_dict['expID'])+'_'+txtPoint+saving_format)
            
            if pointNum==0:
                save_image(PropagImageTmp[ibin], expPathEn[ibin]+'propag/PropagImage_'+str(exp_dict['expID'])+'_'+saving_format)
                save_image(White[ibin], expPathEn[ibin]+'White_'+str(exp_dict['expID'])+'_'+saving_format)
                

    experiment.saveAllParameters(time0,exp_dict)
    
    print("\nfini")
    



        