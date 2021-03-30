#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:37:59 2020

@author: quenot
"""

from xml.dom import minidom
import sys
import datetime
import time
sys.path.append('CodePython/InputOutput/')
sys.path.append('PhaseRetrieval2020')
sys.path.append('PhaseRetrieval2020/InputOutput')
import os
from InputOutput.pagailleIO import saveEdf,openImage
from Experiment import Experiment
import numpy as np

    
     
if __name__ == "__main__":
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
    exp_dict['simulation_type']="Fresnel" #"Fresnel" or "RayT" 

   
    #************************************************************************
    now=datetime.datetime.now()
    exp_dict['expID']=now.strftime("%Y%m%d-%H%M%S") 
    
    SampleImage=[]
    ReferenceImage=[]
    PropagImage=[]
    AbsImage=[]
    SubImage=[]
    Geometry=[]
    for pointNum in range(exp_dict['nbExpPoints']):
        print("\n\nINITIALIZING EXPERIMENT PARAMETERS AND GEOMETRIES")
        experiment=Experiment(exp_dict, pointNum) 
        print("\n\n*************************")
        print("Calculations point",pointNum)
        if exp_dict['simulation_type']=="Fresnel":
            SampleImageTmp, ReferenceImageTmp,PropagImageTmp, White=experiment.computeSampleAndReferenceImages(exp_dict)
        elif exp_dict['simulation_type']=="RayT":
            SampleImageTmp, ReferenceImageTmp,PropagImageTmp, White, Dx, Dy=experiment.computeSampleAndReferenceImagesRT(exp_dict)
        else:
            raise Exception("simulation Type not defined: ", exp_dict['simulation_type'])
        Nbin=len(SampleImageTmp)
        # White=White[exp_dict['margin']:-exp_dict['margin'],exp_dict['margin']:-exp_dict['margin']]
        if pointNum==0:
            PropagImage.append(PropagImageTmp)#[exp_dict['margin']:-exp_dict['margin'],exp_dict['margin']:-exp_dict['margin']])#/White)
        SampleImage.append(SampleImageTmp)#[exp_dict['margin']:-exp_dict['margin'],exp_dict['margin']:-exp_dict['margin']])#/White)
        ReferenceImage.append(ReferenceImageTmp)#[exp_dict['margin']:-exp_dict['margin'],exp_dict['margin']:-exp_dict['margin']])#/White)
        SubImage.append(ReferenceImage[pointNum]-SampleImage[pointNum]/PropagImage[0])
        Geometry.append(experiment.myMembrane.myGeometry[0])

    acquisitions={}
    acquisitions["SampleImage"]=np.asarray(SampleImage)
    acquisitions["ReferenceImage"]=np.asarray(ReferenceImage)
    acquisitions["PropagImage"]=np.asarray(PropagImage)
    acquisitions["SubImage"]=np.asarray(SubImage)
    acquisitions["AbsImage"]=np.asarray(AbsImage)
    acquisitions["White"]=np.asarray(White)

    exp_dict=experiment.createExpDict(exp_dict)
    
    ##    SAVING
    if save: #Creates a text file with all experiment parameters
        print("\n***************************************************")
        print("SAVING\n")
        if exp_dict['simulation_type']=="Fresnel":
            expImagesFilePath=exp_dict['filepath']+'Fresnel_'+str(exp_dict['expID'])+'/'
        if exp_dict['simulation_type']=="RayT":
            expImagesFilePath=exp_dict['filepath']+'RayTracing_'+str(exp_dict['expID'])+'/'
        os.mkdir(expImagesFilePath)
        os.mkdir(expImagesFilePath+'ref/')
        os.mkdir(expImagesFilePath+'sample/')
        os.mkdir(expImagesFilePath+'propag/')
        os.mkdir(expImagesFilePath+'sub/')
        os.mkdir(expImagesFilePath+'membraneThickness/')
#    
        saveEdf(PropagImage[0], expImagesFilePath+'propag/PropagImage_'+str(exp_dict['expID'])+'_'+'.edf')
        saveEdf(White, expImagesFilePath+'White_'+str(exp_dict['expID'])+'_'+'.edf')
#   
        for pointNum in range(exp_dict['nbExpPoints']):
            txtPoint = '%2.2d' % pointNum
            
            saveEdf(Geometry[pointNum], expImagesFilePath+'membraneThickness/'+exp_dict['experimentName']+'_sampling'+str(exp_dict['sampleSampling'])+'_'+str(pointNum)+'.edf')
            saveEdf(SampleImage[pointNum], expImagesFilePath+'sample/sampleImage_'+str(exp_dict['expID'])+'_'+txtPoint+'.edf')
            saveEdf(ReferenceImage[pointNum], expImagesFilePath+'ref/ReferenceImage_'+str(exp_dict['expID'])+'_'+txtPoint+'.edf')
            saveEdf(SubImage[pointNum], expImagesFilePath+'sub/SubImage_'+str(exp_dict['expID'])+'_'+txtPoint+'.edf')
        
        experiment.saveAllParameters(time0,exp_dict)
    
    print("\nfini")


        