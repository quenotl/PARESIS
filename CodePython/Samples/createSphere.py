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

    Sample=np.zeros((dimX,dimY))
    
    myRadius=myRadius/pixelSize
    
    for i in range(dimX):
        for j in range(dimY):
            dist=(dimY/2-i)**2+(dimY/2-j)**2
            if dist<myRadius**2:
                Sample[i,j]=2*np.sqrt(myRadius**2-(dimY/2-j)**2-(dimX/2-i)**2)
            
#    
#    print('\nSample sohere geom')
#    plt.figure()
#    plt.imshow(Sample*pixelSize*1e-6)
#    plt.colorbar()
#    plt.show()
            
    return Sample*pixelSize*1e-6
            

def getText(node):
    return node.childNodes[0].nodeValue