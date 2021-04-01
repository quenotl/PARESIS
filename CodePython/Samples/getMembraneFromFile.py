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
import cv2
import os

def getMembraneFromFile(myMembraneFile,studyDimensions,studyPixelSize, oversamp, numPoint):
    
    fullImagePaths = glob.glob(myMembraneFile + '/*.tif') + glob.glob(myMembraneFile + '/*.tiff') + glob.glob(myMembraneFile + '/*.edf')
    fullImagePaths.sort()
    thickness=(openImage(fullImagePaths[numPoint]))
    Nx,Ny=thickness.shape
    
    
    if studyDimensions[0]!=Nx or studyDimensions[1]!=Ny:
        
        raise ValueError("The membrane you are trying to load does not have the correct dimensions")
               
    thickness = np.asarray(thickness, dtype=np.float)
    
    print('\Membrane Geometry')
    plt.figure()
    plt.imshow(thickness)
    plt.colorbar()
    plt.show()
    
    return (thickness)


def getMembraneSegmentedFromFile(sample,dimX,dimY,pixSize,overSamp, pointNum):
    margin=int(np.ceil(10*sample.myMeanSphereRadius/pixSize))
    margin2=int(np.floor(margin/2))
    if sample.myMembraneFile.split('/')[-1]=='CuSn.txt':
        corrFactor=sample.myMeanSphereRadius/6
        membraneSizeinFiley=int(np.floor(8178))*corrFactor 
        membraneSizeinFilex=int(np.floor(9813))*corrFactor 
    else: 
        raise Exception("Enter segmented membrane size ")
    
    with open(sample.myMembraneFile) as file:
        parameters=json.load(file)
        
    
    membrane=np.zeros((dimX+2*margin,dimY+2*margin))
    Nsphere=len(parameters)
    
    Offsetx=0
    Offsety=0
    offsetCorr=0
    
    for nmem in range(sample.myNbOfLayers):
        if pointNum!=0 or nmem!=0:
            maxOffsetx=int(np.floor((membraneSizeinFilex/pixSize)/2-dimX/2))-margin
            maxOffsety=int(np.floor((membraneSizeinFiley/pixSize)/2-dimY/2))-margin
            Offsetx=np.random.random_integers(-maxOffsetx,maxOffsetx)*pixSize
            Offsety=np.random.random_integers(-maxOffsety,maxOffsety)*pixSize
        
        for i in range(Nsphere):
            radFloat=parameters[i][2]*corrFactor/pixSize
            radInt=int(np.floor(radFloat))+1
            x=int(np.round(parameters[i][1]*corrFactor/pixSize+dimX/2+margin+offsetCorr+Offsety/pixSize))
            y=int(np.round(parameters[i][0]*corrFactor/pixSize+dimY/2+margin+offsetCorr+Offsetx/pixSize))
            xfloat=parameters[i][1]*corrFactor/pixSize+dimX/2+margin+Offsety/pixSize+offsetCorr
            yfloat=parameters[i][0]*corrFactor/pixSize+dimY/2+margin+Offsetx/pixSize+offsetCorr
            if margin2<x<dimX+margin+margin2 and margin2<y<dimY+margin+margin2:
    #            print(x,y,radFloat)
                for ii in range(-radInt,radInt):
                    for jj in range(-radInt,radInt):
                        dist=np.sqrt(((ii+x-xfloat))**2+((jj+y-yfloat))**2)
                        if dist<radFloat:
                            membrane[x+ii,y+jj]+=2*np.sqrt(radFloat**2-dist**2)
                            
    membrane=membrane[margin:margin+dimX,margin:margin+dimY]
    return membrane*pixSize*1e-6



    
if __name__ == "__main__":

    class Membrane():
        def __init__(self):
            self.myMembraneFile='Membranes/CuSn.txt'
            self.myMeanSphereRadius=10
            self.myNbOfLayers=2
            return
    
    # PARAMETERS
    number_of_positions=10
    imageMargins=10
    overSamp=4
    Nx=200 #Detector size
    Ny=200
    dimX=int((Nx+2*imageMargins)*overSamp)
    dimY=int((Ny+2*imageMargins)*overSamp)
    dist_source_sample=0.1 #in m
    dist_sample_detector=0.5 #in m
    detector_pixel_size=24 #in um
    magnification=(dist_sample_detector+dist_source_sample)/dist_source_sample
    pixSize=detector_pixel_size/overSamp/magnification
    storingFolder='Membranes/CuSn_dim'+str(Nx)+'x'+str(Ny)+'_oversamp'+str(overSamp)+'_margin'+str(imageMargins)+'/'
    
    membrane=Membrane()
    
    maxx=(9813/6*membrane.myMeanSphereRadius)
    maxy=(8178/6*membrane.myMeanSphereRadius)
    maxDx=maxx/pixSize
    maxDy=maxy/pixSize

    print("Maximum size of the membrane you can generate with those parameters: %d x %d um" %(maxx ,maxy) )
    print("The one you are trying to create is %d x %d um" %(dimX*pixSize, dimY*pixSize))
    if dimX>maxDx or dimY>maxDy:
        raise Exception("The segmented membrane is too small for the one you are trying to generate. (You can try increasing the sphere mean radius)")
    
    os.mkdir(storingFolder)
    
    for pointNum in range (number_of_positions):
        print("\nNumber of the position currently processed ", pointNum)
        membraneThickness=getMembraneSegmentedFromFile(membrane,dimX,dimY,pixSize,overSamp, pointNum)
        txtPoint = '%2.2d' % pointNum
        saveEdf(membraneThickness, storingFolder+txtPoint+'.edf')
    