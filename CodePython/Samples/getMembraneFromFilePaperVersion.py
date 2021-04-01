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
import cv2
import os

def getMembraneFromFile(myMembraneFile,nbOfLayers,studyDimensions,studyPixelSize, oversamp, numPoint):
    fullImage=(openImage(myMembraneFile))
    Nx,Ny=fullImage.shape
    thickness=np.zeros(studyDimensions)
#    resizeFactor=0.65
#    Nx2=int(np.floor(Nx*resizeFactor))
#    Ny2=int(np.floor(Ny*resizeFactor))
#    fullImage=(cv2.resize(fullImage,(Ny,Nx)))
#    Nx,Ny=fullImage.shape
    
    
    if studyDimensions[0]>Nx or studyDimensions[1]>Ny:
        
        raise ValueError("The membrane you are trying to load is too small for study dimensions")
    
    for iSup in range(nbOfLayers):
        if numPoint!=0 or iSup!=0:
            Offsetx=np.random.random_integers(0,int(np.round(Nx-studyDimensions[0])))#numPoint*stepSize%4
            Offsety=np.random.random_integers(0,int(np.round(Ny-studyDimensions[1])))#int(np.floor(numPoint*stepSize/4))
            if numPoint%2==1 or iSup%2==1:
                fullImage=np.rot90(np.rot90(fullImage))
        else:
            Offsetx=int(np.round((Nx-studyDimensions[0])/2))
            Offsety=int(np.round((Ny-studyDimensions[1])/2))
    #    for i in range(numPoint):
    #        fullImage=np.rot90(fullImage)
                    
        thickness+=fullImage[Offsetx:studyDimensions[0]+Offsetx,Offsety:Offsety+studyDimensions[1]]
    
    thickness = np.asarray(thickness, dtype=np.float)
    
    #thickness2=thickness*studyPixelSize*1e-6*factor*resizeFactor*5
    
    print('\Membrane Geometry')
    plt.figure()
    plt.imshow(thickness)
    plt.colorbar()
    plt.show()
    
    return (thickness)

def getMembraneFromMultipleFiles(myMembraneFolder,nbOfLayers,studyDimensions,studyPixelSize, oversamp, numPoint):
    refImages = glob.glob(myMembraneFolder + '/*.edf')
    
    fullImage=(openImage(refImages))
    Nx,Ny=fullImage.shape
    thickness=np.zeros(studyDimensions)
    if studyDimensions[0]>Nx or studyDimensions[1]>Ny:
        
        raise ValueError("The membrane you are trying to load is too small for study dimensions")
    
    for iSup in range(nbOfLayers):
        Offsetx=int(np.round((Nx-studyDimensions[0])/2))
        Offsety=int(np.round((Ny-studyDimensions[1])/2))
    #    for i in range(numPoint):
    #        fullImage=np.rot90(fullImage)
                    
        thickness+=fullImage[Offsetx:studyDimensions[0]+Offsetx,Offsety:Offsety+studyDimensions[1]]
    
    thickness = np.asarray(thickness, dtype=np.float)
    
    #thickness2=thickness*studyPixelSize*1e-6*factor*resizeFactor*5
    
    print('\Membrane Geometry')
    plt.figure()
    plt.imshow(thickness)
    plt.colorbar()
    plt.show()
    
    return (thickness)
    


def getMembraneSegmentedFromFile(sample,dimX,dimY,pixSize,overSamp, pointNum):
    margin=int(np.ceil(10*sample.myMeanSphereRadius/pixSize))#int(np.floor(500/pixSize))
    margin2=int(np.floor(margin/2))
    if sample.myMembraneFile.split('/')[-1]=='Fe90_3D.txt':
        membraneSizeinFiley=int(np.floor(12244.8*3)) #11313 for Fe90
        membraneSizeinFilex=int(np.floor(11313.6)) #11313 for Fe90
    elif sample.myMembraneFile.split('/')[-1]=='CuSn.txt':
        corrFactor=sample.myMeanSphereRadius/6
        membraneSizeinFiley=int(np.floor(8178))*corrFactor #11313 for Fe90
        membraneSizeinFilex=int(np.floor(9813))*corrFactor #11313 for Fe90
        listOffsetsy=[0,200,400,400,400]
        listOffsetsx=[0,0,0,200,400]
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
        if sample.myMembraneFile.split('/')[-1]=='CuSn.txt':
            Offsetx=listOffsetsx[pointNum]
            Offsety=listOffsetsy[pointNum]
            if nmem!=0:
                print("OFFSET DONE")
                corrFactor=41.5/44.5
                Offsetx+=190*corrFactor
        # if pointNum!=0 or nmem!=0:
        #     maxOffsetx=int(np.floor((membraneSizeinFilex/pixSize)/2-dimX/2))-margin
        #     maxOffsety=int(np.floor((membraneSizeinFiley/pixSize)/2-dimY/2))-margin
        #     Offsetx=np.random.random_integers(-maxOffsetx,maxOffsetx)*pixSize#numPoint*stepSize%4
        #     Offsety=np.random.random_integers(-maxOffsety,maxOffsety)*pixSize#int(np.floor(numPoint*stepSize/4))
        
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
#    filename="/Users/quenot/Documents/Simulations/CodePython_2605_membrane/Samples/Membranes/membraneThickness_poro0.7_rad22.5_sig6"
#    filename='/Users/quenot/Documents/Simulations/CodePython_18012021/Samples/Membranes/Fe90_3D.txt'
#    


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
    # os.mkdir(storingFolder)
    membrane=Membrane()

    
    for pointNum in range (number_of_positions):
        membraneThickness=getMembraneSegmentedFromFile(membrane,dimX,dimY,pixSize,overSamp, pointNum)
        txtPoint = '%2.2d' % pointNum
        saveEdf(membraneThickness, storingFolder+txtPoint+'.edf')
    