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
    exp_dict['experimentName']="SIMAP_FilNylon_mars2021"
    # Output filepath to store the result images
    exp_dict['filepath']='../Results/SIMAP_FilNylon_mars2021/'
    # Define algorithm parameters
    exp_dict['sampleSampling']=4 # MUST BE AN INTEGER
    exp_dict['nbExpPoints']=1 #number of pair of acquisitions (Ir, Is) simulated with different positions of the membrane
    exp_dict['margin']=10 #with Fresnel there might be an aliasing issue so we need to extend very slightly the image for calculations
    save=True
    exp_dict['simulation_type']="RayT" #"Fresnel" or "RayT" 

   
    #************************************************************************
    now=datetime.datetime.now()
    exp_dict['expID']=now.strftime("%Y%m%d-%H%M%S") #define experiment ID
    
    SampleImage=[]
    ReferenceImage=[]
    PropagImage=[]
    AbsImage=[]
    SubImage=[]
    Geometry=[]
    
    experiment=Experiment(exp_dict) 
    
    for pointNum in range(exp_dict['nbExpPoints']):
        experiment.myMembrane.myGeometry=[]
        experiment.myMembrane.getMyGeometry(experiment.studyDimensions,experiment.myMembrane.membranePixelSize,experiment.sampling, pointNum, exp_dict['nbExpPoints'])
        print("\n\nINITIALIZING EXPERIMENT PARAMETERS AND GEOMETRIES")
        print("\n\n*************************")
        print("Calculations point",pointNum)
        if exp_dict['simulation_type']=="Fresnel":
            SampleImageTmp, ReferenceImageTmp,PropagImageTmp, White=experiment.computeSampleAndReferenceImages(pointNum)
        elif exp_dict['simulation_type']=="RayT":
            SampleImageTmp, ReferenceImageTmp,PropagImageTmp, White, Dx, Dy=experiment.computeSampleAndReferenceImagesRT(pointNum)
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
                os.mkdir(expPathEn[ibin]+'ref/')
                os.mkdir(expPathEn[ibin]+'sample/')
                os.mkdir(expPathEn[ibin]+'propag/')
            
        txtPoint = '%2.2d' % pointNum
        saveEdf(experiment.myMembrane.myGeometry[0], expImagesFilePath+'membraneThickness/'+exp_dict['experimentName']+'_sampling'+str(exp_dict['sampleSampling'])+'_'+str(pointNum)+'.edf')
            
        for ibin in range(Nbin):
            saveEdf(SampleImageTmp[ibin], expPathEn[ibin]+'sample/sampleImage_'+str(exp_dict['expID'])+'_'+txtPoint+'.edf')
            saveEdf(ReferenceImageTmp[ibin], expPathEn[ibin]+'ref/ReferenceImage_'+str(exp_dict['expID'])+'_'+txtPoint+'.edf')
            
            if pointNum==0:
                saveEdf(PropagImageTmp[ibin], expPathEn[ibin]+'propag/PropagImage_'+str(exp_dict['expID'])+'_'+'.edf')
                saveEdf(White[ibin], expPathEn[ibin]+'White_'+str(exp_dict['expID'])+'_'+'.edf')

    experiment.saveAllParameters(time0,exp_dict)
    
    print("\nfini")
    



        