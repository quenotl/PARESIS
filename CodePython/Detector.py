#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:40:55 2020

@author: quenot
"""

from xml.dom import minidom
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import convolve
from matplotlib import pyplot as plt
import time
from numba import jit
import xlrd

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
        self.myEnergyLimit=200 #No longer useful?
        self.myBinsThersholds=[] #keX
        self.myScintillatorMaterial=None
        self.myScintillatorThickness=0. #um
        self.beta=[]
        
        
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
                    if node.localName=="myBinsThersholds":
                        myBinsThersholdsTmp=self.getText(currentDetector.getElementsByTagName("myBinsThersholds")[0])
                        myBinsThersholdsTmp=list(myBinsThersholdsTmp.split(","))
                        self.myBinsThersholds=[float(ele) for ele in myBinsThersholdsTmp]
                    if node.localName=="myScintillatorMaterial":
                        self.myScintillatorMaterial=(self.getText(currentDetector.getElementsByTagName("myScintillatorMaterial")[0]))
                        self.myScintillatorThickness=float(self.getText(currentDetector.getElementsByTagName("myScintillatorThickness")[0]))
                return
            
        raise ValueError("detector not found in xml file")
            
            
    def detection(self,incidentWave,effectiveSourceSize,samplingFactor):
        
        """
        Adds source and PSF blurrings, resamples to detector pixel size and add shot noise

        Args:
            incidentWave (2d numpy array): Intensity arriving to the detector.
            effectiveSourceSize (float): source projected FWHM.

        Returns:
            detectedImage (2d numpy array): detected image.

        """
        sampledpixeledsize=self.myPixelSize/samplingFactor
        weights=self.source_geometry(5e-4,50e-6,sampledpixeledsize)
        if effectiveSourceSize!=0:
            sigmaSource=effectiveSourceSize/2.355 #from FWHM to std dev
            # incidentWave=gaussian_filter(incidentWave, sigmaSource,mode='wrap')
            incidentWave=convolve(incidentWave, weights, mode='nearest')
        intensityBeforeDetection=resize(incidentWave, self.myDimensions[0],self.myDimensions[1])
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
    
    
    def getBeta(self, sourceSpectrum):
        print("Materials :", self.myScintillatorMaterial)
        pathTablesDeltaBeta ='Samples/DeltaBeta/TablesDeltaBeta.xls'
        for sh in xlrd.open_workbook(pathTablesDeltaBeta).sheets():
            for col in range(sh.ncols):
                row=0
                myCell = sh.cell(row, col)
#                    print("\n\n Sample materials : ",self.myMaterials[imat])
                if myCell.value == self.myScintillatorMaterial:
                    row=row+3
                    for energy,_ in sourceSpectrum:
                        currentCellValue=sh.cell(row,col).value
                        if energy*1000<currentCellValue:
                            self.beta.append((energy,1))
                            print("No delta beta values under",currentCellValue, "eV")
                            continue
                        nextCellValue=sh.cell(row+1,col).value
                        while nextCellValue<energy*1e3: #find in which interval the energy is in the delta beta tables
                            row+=1
                            currentCellValue=sh.cell(row,col).value
                            nextCellValue=sh.cell(row+1,col).value
                        #Linear interpollation between those values
                        step=nextCellValue-currentCellValue
                        currentCellBeta=float(sh.cell(row,col+2).value)
                        nextCellBeta=float(sh.cell(row+1,col+2).value)
                        betaInterp=abs(nextCellValue-energy*1e3)/step*currentCellBeta+abs(currentCellValue-energy*1e3)/step*nextCellBeta
                        self.beta.append((energy,betaInterp))
                        
                    return
            raise ValueError("The scintillator material has not been found in delta beta tables")
    
    def source_geometry(self, distance, sigma,pixel_size):
    
            pixel_size=pixel_size*1e-6 
            x = 5*sigma/pixel_size
            y = distance/pixel_size+5*sigma/pixel_size
            y0, x0=y/2-distance/pixel_size/2, x/2
            y1, x1=y/2+distance/pixel_size/2, x/2
        
            x=np.arange(x)
            y=np.arange(y)
            sigma /= pixel_size
    
            gx = np.exp(-(x-x0)**2/(2*sigma**2))
            gx1 = np.exp(-(x-x1)**2/(2*sigma**2))

            gy = np.exp(-(y-y0)**2/(2*sigma**2))
            gy1 = np.exp(-(y-y1)**2/(2*sigma**2))

            g = np.outer(gx, gy)
            g1 = np.outer(gx1, gy1)
            g2=(g+g1)/2
            g=g2/np.sum(g2)
            plt.figure()
            plt.imshow(g2)
            plt.colorbar()
            plt.show()
            return g2


      
    
@jit(nopython=True)
def resize(imageToResize,sizeX, sizeY):
    Nx, Ny=imageToResize.shape
    if Nx==sizeX and Ny==sizeY:
        return imageToResize
    
    resizedImage=np.ones((sizeX,sizeY))
    sampFactor=int(Nx/sizeX)
    
    for x0 in range(sizeX):
        for y0 in range(sizeY):
            resizedImage[x0,y0]=np.sum(imageToResize[int(np.floor(x0*sampFactor)):int(np.floor(x0*sampFactor+sampFactor)),int(np.floor(y0*sampFactor)):int(np.floor(y0*sampFactor+sampFactor))])
            
    return resizedImage


    



