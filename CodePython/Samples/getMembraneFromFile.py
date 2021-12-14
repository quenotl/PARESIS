#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:47:44 2020

@author: quenot
"""

import numpy as np
from InputOutput.pagailleIO import saveEdf,openImage
from matplotlib import pyplot as plt
import json
from InputOutput.pagailleIO import saveEdf,openImage
import glob
import os
import time

def getMembraneFromFile(myMembraneFile,studyDimensions,studyPixelSize, oversamp, numPoint):
    
    fullImagePaths = glob.glob(myMembraneFile + '/*.tif') + glob.glob(myMembraneFile + '/*.tiff') + glob.glob(myMembraneFile + '/*.edf')
    fullImagePaths.sort()
    thickness=(openImage(fullImagePaths[numPoint]))
    if not studyDimensions == thickness.shape:
        raise ValueError("The membrane you are trying to load does not have the correct dimensions")
               
    thickness = np.asarray(thickness, dtype=np.float)
    
    print('\Membrane Geometry')
    plt.figure()
    plt.imshow(thickness)
    plt.colorbar()
    plt.show()
    
    return (thickness)


def getMembraneSegmentedFromFile(sample,dimX,dimY,pixSize,overSamp, pointNum):
    start = time.time()
    margin=int(-(10*sample.myMeanSphereRadius/pixSize//-1)) #in pix Ceiling operator
    margin2=int(margin/2)
    
    if sample.myMembraneFile.split('/')[-1]=='CuSn.txt':
        corrFactor=sample.myMeanSphereRadius/4.57 #to adapt the scale to the mean sphere radius desired
        membraneSizeinFilex=8102*corrFactor+sample.myMeanSphereRadius #??? int(np.floor(8102))
        membraneSizeinFiley=9740*corrFactor+sample.myMeanSphereRadius #??? int(np.floor(9740))
    else: 
        raise Exception("Enter segmented membrane size ")
    
    with open(sample.myMembraneFile) as file:
        parameters=json.load(file)
    parameters=np.asarray(parameters) #change list to numpy array 2D [nsphere, values] 
    parameters=parameters*corrFactor  #to adapt the scale to the mean sphere radius desired
    parameters[:,1]+=membraneSizeinFilex/2 #place origin on top left instead of middle
    parameters[:,0]+=membraneSizeinFiley/2
    
    parameters0=np.copy(parameters)
    membraneSizeinFilex0=membraneSizeinFilex
    membraneSizeinFiley0=membraneSizeinFiley
    
    DiffFileExpX=membraneSizeinFilex/pixSize-dimX
    while DiffFileExpX<0: #as long as the membrane in file is too small, expand it
        print("segmented membrane too small: proceeding with stitching along x")
        parametersStitchX=np.copy(parameters0)
        parametersStitchX[:,1]+=membraneSizeinFilex
        parameters=np.concatenate((parameters, parametersStitchX), axis=0)
        membraneSizeinFilex+=membraneSizeinFilex0
        DiffFileExpX=membraneSizeinFilex/pixSize-dimX
        
    parameters0=np.copy(parameters)
    DiffFileExpY=membraneSizeinFiley/pixSize-dimY
    while DiffFileExpY<0:
        print("segmented membrane too small: proceeding with stitching along y")
        parametersStitchY=np.copy(parameters0)
        parametersStitchY[:,0]+=membraneSizeinFiley
        parameters=np.concatenate((parameters, parametersStitchY), axis=0)
        membraneSizeinFiley+=membraneSizeinFiley0
        DiffFileExpY=membraneSizeinFiley/pixSize-dimY
        
    membrane=np.zeros((dimX+2*margin,dimY+2*margin))
    Nsphere,_=parameters.shape
    
    Offsetx=0
    Offsety=0
    offsetCorr=0
    
    for _ in range(sample.myNbOfLayers):
        # if pointNum!=0 or nmem!=0:
        maxOffsetx=membraneSizeinFilex/pixSize-dimX
        maxOffsety=membraneSizeinFiley/pixSize-dimY
        Offsetx=np.random.randint(margin2,maxOffsetx-margin2)
        Offsety=np.random.randint(margin2,maxOffsety-margin2)
        
        for i in range(Nsphere):
            radFloat=parameters[i,2]/pixSize
            radInt=int(-(radFloat//-1))
            xfloat=parameters[i,1]/pixSize-Offsetx+offsetCorr
            yfloat=parameters[i,0]/pixSize-Offsety+offsetCorr
            x=round(xfloat)
            y=round(yfloat)
            if margin2<x<dimX+margin+margin2 and margin2<y<dimY+margin+margin2:
                for ii in range(-radInt,radInt):
                    for jj in range(-radInt,radInt):
                        dist=np.sqrt(((ii+x-xfloat))**2+((jj+y-yfloat))**2)
                        if dist<radFloat:
                            membrane[x+ii,y+jj]+=2*np.sqrt(radFloat**2-dist**2)

    membrane=membrane[margin:-margin,margin:-margin]
    stop = time.time()-start
    print(f'IT took {stop} seconds')
    return membrane*pixSize*1e-6


if __name__ == "__main__":

    class Membrane():
        def __init__(self):
            self.myMembraneFile='Membranes/CuSn.txt'
            self.myMeanSphereRadius=30
            self.myNbOfLayers=2
            return
    
    # PARAMETERS
    number_of_positions=2
    imageMargins=10
    overSamp=2
    Nx=500 #Detector size
    Ny=1500
    dimX=int((Nx+2*imageMargins)*overSamp)
    dimY=int((Ny+2*imageMargins)*overSamp)
    dist_source_membrane=144 #in m
    dist_membrane_detector=5.2 #in m
    detector_pixel_size=24 #in um
    magnification=(dist_membrane_detector+dist_source_membrane)/dist_source_membrane
    pixSize=detector_pixel_size/overSamp/magnification
    
    membrane=Membrane()
    storingFolder='Membranes/CuSn_dim'+str(Nx)+'x'+str(Ny)+'_oversamp'+str(overSamp)+'_margin'+str(imageMargins)+'_radius'+str(membrane.myMeanSphereRadius)+'/'
    
    if not os.path.exists(storingFolder):
        os.mkdir(storingFolder)
    
    for pointNum in range (number_of_positions):
        print("\nNumber of the position currently processed ", pointNum)
        membraneThickness=getMembraneSegmentedFromFile(membrane,dimX,dimY,pixSize,overSamp, pointNum) 
        plt.figure()
        plt.imshow(membraneThickness)
        plt.colorbar()
        plt.show()
        txtPoint = '%2.2d' % pointNum
        saveEdf(membraneThickness, storingFolder+txtPoint+'.edf')
    