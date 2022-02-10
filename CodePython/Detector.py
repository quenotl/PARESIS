#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:40:55 2020

@author: quenot
"""

from xml.dom import minidom
import numpy as np
from scipy.ndimage.filters import gaussian_filter, median_filter
from scipy.signal import fftconvolve
from matplotlib import pyplot as plt
from getk import getk
import time
from numba import jit
import xlrd

class Detector:
    def __init__(self, exp_dict):
        self.xmlDetectorFileName="xmlFiles/Detectors.xml"
        self.xmldocDetector = minidom.parse(self.xmlDetectorFileName)
        self.myName=""
        self.det_param={}
        self.det_param["myDimensions"]=(0,0)
        self.det_param["myPixelSize"]=0. #en um
        self.det_param["myPSF"]=0. #en pixel
        self.det_param["myBinsThersholds"]=[] #keV
        self.det_param["myScintillatorMaterial"]=None
        self.det_param["myScintillatorThickness"]=0. #um
        self.det_param["photonCounting"]=True
        self.mySpectralEfficiency=[]
        self.beta=[]
        
        # parameters units
        self.det_param["myDimensions_unit"]="pixels"
        self.det_param["myPixelSize_unit"]="um"
        self.det_param["myPSF_unit"]="pixels"
        self.det_param["myBinsThersholds_unit"]="keV" #keV
        self.det_param["myScintillatorThickness_unit"]="um" #um
        
        
        
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
                self.det_param["myDimensions"]=self.getMyDimensions(currentDetector)
                self.det_param["myPixelSize"]=float(self.getText(currentDetector.getElementsByTagName("myPixelSize")[0]))
                self.det_param["myPSF"]=float(self.getText(currentDetector.getElementsByTagName("myPSF")[0]))
                
                for node in currentDetector.childNodes:
                    if node.localName=="myEnergyLimit":
                        self.myEnergyLimit=float(self.getText(currentDetector.getElementsByTagName("myEnergyLimit")[0]))
                    if node.localName=="photonCounting":
                        self.det_param["photonCounting"]=bool(self.getText(currentDetector.getElementsByTagName("photonCounting")[0]))
                    if node.localName=="myBinsThersholds":
                        myBinsThersholdsTmp=self.getText(currentDetector.getElementsByTagName("myBinsThersholds")[0])
                        myBinsThersholdsTmp=list(myBinsThersholdsTmp.split(","))
                        self.det_param["myBinsThersholds"]=[float(ele) for ele in myBinsThersholdsTmp]
                    if node.localName=="myScintillatorMaterial":
                        self.det_param["myScintillatorMaterial"]=(self.getText(currentDetector.getElementsByTagName("myScintillatorMaterial")[0]))
                        self.det_param["myScintillatorThickness"]=float(self.getText(currentDetector.getElementsByTagName("myScintillatorThickness")[0]))
                return
            
        raise ValueError("detector not found in xml file")
            
            
    def detection(self,incidentWave,effectiveSourceSize, exp_param):
        """
        Adds source and PSF blurrings, resamples to detector pixel size and add shot noise

        Args:
            incidentWave (2d numpy array): Intensity arriving to the detector.
            effectiveSourceSize (float): source projected FWHM.

        Returns:
            detectedImage (2d numpy array): detected image.

        """
        #Add margin to avoid edges effect when convolving with PSF and source blurring
        margins=15
        incidentWave=np.pad(incidentWave, margins*exp_param['overSampling'], mode='reflect')

        # convolve with source blurring
        if effectiveSourceSize!=0:
            sigmaSource=effectiveSourceSize/2.355 #from FWHM to std dev
            gaussian=create_gaussian_shape(sigmaSource)
            incidentWave=fftconvolve(incidentWave, gaussian,mode='same')
            
        # incidentWave[incidentWave<0]=0
        # reoverSampling at detector pixel size
        intensityBeforeDetection=resize(incidentWave, self.det_param["myDimensions"][0]+margins*2,self.det_param["myDimensions"][1]+margins*2)
        
        # blurring with detector PSF
        if self.det_param["myPSF"]!=0:
            gaussian=create_gaussian_shape(self.det_param["myPSF"])
            detectedImage=fftconvolve(intensityBeforeDetection, gaussian,mode='same')
        else:
            detectedImage=intensityBeforeDetection
            
        # adding shot noise
        seed = int(np.floor(time.time()*100%(2**32-1))) # no idea why this one is so complicated ... 
        rs = np.random.RandomState(seed)
        detectedImage = rs.poisson((detectedImage))
        
        #Remove margins
        detectedImage=detectedImage[margins:self.det_param["myDimensions"][0]+margins,margins:self.det_param["myDimensions"][1]+margins]
        return detectedImage
        
    
    def getText(self,node):
        return node.childNodes[0].nodeValue
    
    def getMyDimensions(self,node):
        dimX=int(self.getText(node.getElementsByTagName("dimX")[0]))
        dimY=int(self.getText(node.getElementsByTagName("dimY")[0]))
        return np.array([dimX ,dimY])
    
    
    def getBeta(self, sourceSpectrum):
        # print("Materials :", self.det_param["myScintillatorMaterial"])
        pathTablesDeltaBeta ='Samples/DeltaBeta/TablesDeltaBeta.xls'
        for sh in xlrd.open_workbook(pathTablesDeltaBeta).sheets():
            for col in range(sh.ncols):
                row=0
                myCell = sh.cell(row, col)
#                    print("\n\n Sample materials : ",self.myMaterials[imat])
                if myCell.value == self.det_param["myScintillatorMaterial"]:
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
            
            
    def getSpectralEfficiency(self):
        i=0
        plotEff=[]
        plotEn=[]
        for energyData, betaEn in self.beta:
            k=getk(energyData*1000)
            currEfficiency=1-np.exp(-2*k*self.det_param["myScintillatorThickness"]*1e-6*betaEn)
            self.mySpectralEfficiency.append((energyData, currEfficiency))
            plotEff.append(currEfficiency)
            plotEn.append(energyData)
            i+=1
        if len(plotEff)>1:
            plt.figure()
            plt.plot(plotEn,plotEff)
            plt.xlabel('Energy (keV)')
            plt.ylabel('Attenuation')
            plt.title("Detector scintillator attenuation power")
            plt.show()
        else:
            print(f'Scintillator attenuation at {plotEn[0]} keV: {plotEff[0]}')
                
    
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


def create_gaussian_shape(sigma):
    dim=round(sigma*3)*2+1
    
    Qx, Qy = np.meshgrid((np.arange(0, dim) - np.floor(dim / 2) ), (np.arange(0, dim) - np.floor(dim / 2)) ) #frequency ranges of the images in fqcy space

    #sigmaY = sig_scale

    g = np.exp(-(((Qx)**2) / 2. / sigma**2 + ((Qy)**2) / 2. / sigma**2))

    return g/np.sum(g)