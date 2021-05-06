#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:40:55 2020

@author: quenot
"""

from xml.dom import minidom
import numpy as np
from scipy.ndimage.filters import gaussian_filter, median_filter
from matplotlib import pyplot as plt
import time
from numba import jit

class Detector:
    def __init__(self, exp_dict):
        self.xmlDetectorFileName="xmlFiles/Detectors.xml"
        self.xmldocDetector = minidom.parse(self.xmlDetectorFileName)
        self.myName=""
        self.myDimensions=(0,0)
        self.myPixelSize=0. #en um
        self.myPSF=0. #en pixel
        self.myEfficiencyLimit=100. #en kev
        self.margins=exp_dict['margin']
        self.myEnergyLimit=200
        
        
    def defineCorrectValuesDetector(self):
        """
        Define detector parameters from the xml file

        Raises:
            ValueError: detector not found in xml file.

        Returns:
            None.

        """
        for currentDetector in self.xmldocDetector.documentElement.getElementsByTagName("detector"):
            correctDetector = self.getText(currentDetector.getElementsByTagName("name")[0])
            if correctDetector == self.myName:
                self.myDimensions=self.getMyDimensions(currentDetector)
                self.myPixelSize=float(self.getText(currentDetector.getElementsByTagName("myPixelSize")[0]))
                self.myPSF=float(self.getText(currentDetector.getElementsByTagName("myPSF")[0]))
                
                for node in currentDetector.childNodes:
                    if node.localName=="myEnergyLimit":
                        self.myEnergyLimit=float(self.getText(currentDetector.getElementsByTagName("myEnergyLimit")[0]))
            
                return
            
        raise ValueError("detector not found in xml file")
            
            
    def detection(self,incidentWave,effectiveSourceSize):
        """
        Adds source and PSF blurrings, resamples to detector pixel size and add shot noise

        Args:
            incidentWave (2d numpy array): Intensity arriving to the detector.
            effectiveSourceSize (float): source projected FWHM.

        Returns:
            detectedImage (2d numpy array): detected image.

        """
        if effectiveSourceSize!=0:
            sigmaSource=effectiveSourceSize/2.355 #from FWHM to std dev
            incidentWave=gaussian_filter(incidentWave, sigmaSource,mode='wrap')
        intensityBeforeDetection=self.resize(incidentWave, self.myDimensions[0],self.myDimensions[1])
        seed       = int(np.floor(time.time()*100%(2**32-1)))
        rs         = np.random.RandomState(seed)
        if self.myPSF!=0:
            detectedImage=gaussian_filter(intensityBeforeDetection, self.myPSF,mode='wrap')
        else:
            detectedImage=intensityBeforeDetection
        
        detectedImage = rs.poisson((detectedImage))

        detectedImage=detectedImage[self.margins:self.myDimensions[0]-self.margins,self.margins:self.myDimensions[1]-self.margins]
        
        return detectedImage
        
    
    def getText(self,node):
        return node.childNodes[0].nodeValue
    
    def getMyDimensions(self,node):
        dimX=int(self.getText(node.getElementsByTagName("dimX")[0]))+self.margins*2
        dimY=int(self.getText(node.getElementsByTagName("dimY")[0]))+self.margins*2
        return np.array([dimX ,dimY])
    
    @jit(nopython=True)
    def resize(self,imageToResize,sizeX, sizeY):
        Nx, Ny=imageToResize.shape
        if Nx==sizeX and Ny==sizeY:
            return imageToResize
        
        resizedImage=np.ones((sizeX,sizeY))
        sampFactor=int(Nx/sizeX)
        
        for x0 in range(sizeX):
            for y0 in range(sizeY):
                resizedImage[x0,y0]=np.sum(imageToResize[int(np.floor(x0*sampFactor)):int(np.floor(x0*sampFactor+sampFactor)),int(np.floor(y0*sampFactor)):int(np.floor(y0*sampFactor+sampFactor))])
                
        return resizedImage