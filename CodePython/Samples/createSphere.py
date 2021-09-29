#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 17:37:42 2020

@author: quenot
"""
import numpy as np
from xml.dom import minidom
import random
from InputOutput.pagailleIO import openImage, saveEdf, openSeq
from matplotlib import pyplot as plt

def CreateSampleSphere(myName, dimX, dimY, pixelSize):
    print("Je cree l'echantillon!!")

    xmlSampleFileName="Samples/Samples.xml"
    xmldocSample = minidom.parse(xmlSampleFileName)
    
    for currentSample in xmldocSample.documentElement.getElementsByTagName("sample"):
        correctSample = getText(currentSample.getElementsByTagName("name")[0])
        
        if correctSample == myName:
            
            myRadius=float(getText(currentSample.getElementsByTagName("myRadius")[0]))

    Sample=np.zeros((1,dimX,dimY))
    
    myRadius=myRadius/pixelSize
    
    for i in range(dimX):
        for j in range(dimY):
            dist=(dimY/2-i)**2+(dimY/2-j)**2
            if dist<myRadius**2:
                Sample[0,i,j]=2*np.sqrt(myRadius**2-(dimY/2-j)**2-(dimX/2-i)**2)
            
#    
#    print('\nSample sohere geom')
#    plt.figure()
#    plt.imshow(Sample*pixelSize*1e-6)
#    plt.colorbar()
#    plt.show()
            
    return Sample*pixelSize*1e-6


def CreateSampleSpheresInTube(myName, dimX, dimY, pixelSize):
    print("Je cree l'echantillon!!")
    print("Pixel size", pixelSize)
            
    myRadius=500 #um

    Sample=np.zeros((3,dimX,dimY))
    
    myRadius=myRadius/pixelSize
    patchSize2=int(np.ceil(myRadius))
    patchSize=patchSize2*2
    spherePatch=np.zeros((patchSize,patchSize))
    
    for i in range(patchSize):
        for j in range(patchSize):
            dist=(patchSize/2-i)**2+(patchSize/2-j)**2
            if dist<myRadius**2:
                spherePatch[i,j]=2*np.sqrt(myRadius**2-(patchSize/2-j)**2-(patchSize/2-i)**2)
            
    posX=int(np.round(dimX/2))
    posYmuscle=int(np.round(1500/pixelSize))
    posYcart=int(np.round(3500/pixelSize))
    
    
    print('\nSample sohere geom')
    plt.figure()
    plt.imshow(spherePatch*pixelSize*1e-6)
    plt.colorbar()
    plt.show()

            
    myRadius=1500 #um
    WireAngle=0
    
    Tube=np.zeros((dimX,dimY))
    Nxp=2*dimX
    Nyp=2*dimY
    diffx=int((Nxp-dimX)/2)
    diffy=int((Nyp-dimY)/2)

    myRadius=myRadius/pixelSize
    
    for j in range(dimX):
        if (abs(dimX/2-j)<myRadius):
            Tube[j,:]=2*np.sqrt(myRadius**2-(dimX/2-j)**2)
    
    print("Tube")
    plt.figure()
    plt.imshow(Tube)
    plt.colorbar()
    plt.show()
    
    Sample[0,posX-patchSize2:posX+patchSize2,posYmuscle-patchSize2:posYmuscle+patchSize2]=spherePatch
    Sample[1,posX-patchSize2:posX+patchSize2,posYcart-patchSize2:posYcart+patchSize2]=spherePatch
    Sample[2]=Tube-Sample[0]-Sample[1]
        
    Sample[0]=Sample[0].transpose()
    Sample[1]=Sample[1].transpose()
    Sample[2]=Sample[2].transpose()
            
    return Sample*pixelSize*1e-6
            

def getText(node):
    return node.childNodes[0].nodeValue