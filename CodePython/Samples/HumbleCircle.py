# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 11:55:39 2021

@author: Robin
"""

import numpy as np
from xml.dom import minidom
import random
from InputOutput.pagailleIO import openImage, saveEdf, openSeq
from matplotlib import pyplot as plt

def CreateSampleSphereMultiple(dimX, dimY, pixelSize):
    
    Nspheres=4
    spheresRadius=[2500, 5000, 7500, 10000]
    
    imageCenterHorizontal=int(np.round(dimY/2))
   
    Sample=np.zeros((1,dimX,dimY))

    for ns in range(Nspheres):
        spheresCenters = int(np.round(dimX/(Nspheres*2))) + ns * int(np.round(dimX / Nspheres))
        
        myRadius=(spheresRadius[ns]/pixelSize) #/ 7
        print("Radius en px : ", myRadius)
        
        centerx=spheresCenters
        for i in range(dimX):
            for j in range(dimY):
                dist = (centerx-i)**2+(imageCenterHorizontal-j)**2
                if dist<myRadius**2:
                    Sample[0,i,j]=2*np.sqrt(myRadius**2-(imageCenterHorizontal-j)**2-(centerx-i)**2)

    print('\nMutiple Spheres sample')
    plt.figure()
    plt.imshow(Sample[0]*pixelSize*1e-3)
    plt.colorbar()
    plt.show()

    return Sample*pixelSize*1e-6


def getText(node):
    return node.childNodes[0].nodeValue

if __name__ == "__main__":

    CreateSampleSphereMultiple(5680, 1480, 1.76)