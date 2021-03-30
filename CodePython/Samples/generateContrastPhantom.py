#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 15:05:44 2020

@author: quenot
"""
import numpy as np
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2
from numpy.fft import ifft2 as ifft2
from matplotlib import pyplot as plt
from skimage.transform import radon, rescale
import json
import os
from InputOutput.pagailleIO import saveEdf,openImage
import glob

def generateContrastPhantom(dimX, dimY, pixsize, angle):
    pixsize=pixsize/1000 #m to mm
    geometry=[]
    Slice=np.zeros(((13,dimY,dimY)))
    lines=np.zeros(dimY)
    projTot=np.zeros((dimX, dimY))
    smallTubesRadius=2 #mm
    supportRadius=15
    TubesCenters=[[22,7],[16.4,4.7],[10,7],[6,12.5],[6,19.5],[10,25],[16.4,27],[22,25],[21,16],[11,16],[16,11],[16,21]]
    Nmat=12        
    TubesCenters=np.asarray(TubesCenters, dtype=np.float64) #mm

    
    sliceTotMat=np.zeros((dimY, dimY))
    
    smallTubesRadiuspix=smallTubesRadius/pixsize #pix
    smallTubesRadiuspixInt=int(np.floor(smallTubesRadiuspix)+2) #pixEntier
    supportRadiuspix=supportRadius/pixsize
    supportRadiuspixInt=int(np.floor(supportRadiuspix)+2)
    
    origin=dimY/2-16/pixsize #pix
    
    if origin<0:
        raise ValueError("Image too small for the contrast phantom size")
    
    TubesCenters= TubesCenters/pixsize+origin #pix

    for imat in range(Nmat):
        sliceMat=np.zeros((dimY, dimY))
        centeri=TubesCenters[imat,0] #pix
        centerj=TubesCenters[imat,1]
        centeriInt=int(np.round(centeri)) #pixInt
        centerjInt=int(np.round(centerj))
        for ii in range(centeriInt-smallTubesRadiuspixInt,centeriInt+smallTubesRadiuspixInt):
            for jj in range(centerjInt-smallTubesRadiuspixInt,centerjInt+smallTubesRadiuspixInt):
                dist=((ii-centeri))**2+((jj-centerj))**2
                if dist<smallTubesRadiuspix**2:
                    sliceMat[ii,jj]=1
                    
        Slice[imat]=sliceMat
        sliceTotMat+=sliceMat*(imat+1)
        
    
    sliceMat=np.zeros((dimY, dimY))
    centeri=dimY/2 #pix
    centerj=dimY/2
    centeriInt=int(np.round(centeri)) #pixInt
    centerjInt=int(np.round(centerj))
    for ii in range(centeriInt-supportRadiuspixInt,int(np.round(27/pixsize+origin))):
        for jj in range(centerjInt-supportRadiuspixInt,centerjInt+supportRadiuspixInt):
            if sliceTotMat[ii,jj]==0:
                dist=((ii-centeri))**2+((jj-centerj))**2
                if dist<supportRadiuspix**2:
                    sliceMat[ii,jj]=1
                
    Slice[12]=sliceMat
    sliceTotMat+=sliceMat*(13)
    
    print('Slice Geometry')
    plt.figure()
    plt.imshow(sliceTotMat)
    plt.colorbar()
    plt.show()
    
    for imat in range(Nmat+1):
        
        projMat=np.ones((dimX,dimY))
        Line=radon(Slice[imat],[angle])
        Nx,Ny=Line.shape
        Line=np.transpose(Line[int(np.round(Nx/2-dimY/2)):int(np.round(Nx/2+dimY/2))])
        
        if imat<12:
            projMat=projMat*Line
        if imat==12:
            projMat[int(np.floor(dimX/2)):dimX]=projMat[int(np.floor(dimX/2)):dimX]*Line
            projMat[0:int(np.floor(dimX/2))]=0    
        
        projTot+=projMat*(imat+1)*pixsize*1e-3
        geometry.append(projMat*pixsize*1e-3)
    
    
    print('Proj Geometry')
    plt.figure()
    plt.imshow(geometry[12])
    plt.colorbar()
    plt.show()
    
    return geometry
    
def openContrastPhantom(myGeometryFolder,dimX, dimY, pixsize,oversamp, angle):
    pixsizeStr='%4.1f'%(pixsize*oversamp)
    filepath="Samples/"+myGeometryFolder
    filepaths=glob.glob(filepath+"/*")
    filepaths.sort()
    print(filepath)
    geometry=[]
    if filepaths!=[]:
        for i in range(13):
            geometry.append(openImage(filepaths[i]))
    else:
         raise Exception("The contrast phantom geometry you are trying to load does not exist or is incorrectly named", filepath)   
    
    return geometry
        
        
    
if __name__ == "__main__":
    
    oversamp=1
    Ny=1500*oversamp
    Nx=500*oversamp
    angle =90
    pixelSize=21.4/oversamp #um
    geometry2=generateContrastPhantom(Nx, Ny, pixelSize, angle)
    dim=geometry2[0].shape
#    geometry=geometry.tolist()

    filepath="ContrastPhantom/CP_"+str(pixelSize*oversamp).replace(".","p")+"_angle"+str(angle)
    os.mkdir(filepath)
    for i in range(13):
        saveEdf(geometry2[i], filepath+"/mat"+str(i)+".edf")
