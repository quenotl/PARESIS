#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:37:59 2020

@author: quenot
"""

import numpy as np
from xml.dom import minidom
import time
from Source import Source
from Detector import Detector
from Sample import AnalyticalSample
from numpy.fft import fftshift
from numpy.fft import ifftshift
from numpy.fft import fft2
from numpy.fft import ifft2
from matplotlib import pyplot as plt
from numpy import pi as pi

class Experiment:
    def __init__(self, exp_dict, pointNum):#expName, pointNum, sampling):
        """Experiment class constructor.
            
        Args:
            exp_dict (dict): dictionnary of simulation algorithm parameters
        """
        self.xmlExperimentFileName="xmlFiles/Experiment.xml"
        self.xmldoc = minidom.parse(self.xmlExperimentFileName)
        self.name=exp_dict['experimentName']
        self.distSourceToMembrane=0
        self.distMembraneToObject=0
        self.distObjectToDetector=0
        self.mySampleType=""
        self.sampling=exp_dict['sampleSampling']
        self.studyPixelSize=0.
        self.studyDimensions=(0.,0.)
        self.mySampleofInterest=None
        self.myDetector=None
        self.mySource=None
        self.myMembrane=None
        self.myPlaque=None
        self.myAirVolume=None

        self.meanShotCount=0
        self.meanEnergy=0
        self.nbPoints=pointNum
        self.Dxreal=[]
        self.Dyreal=[]
        self.imageSampleBeforeDetection=[]
        self.imageReferenceBeforeDetection=[]
        self.imagePropagBeforeDetection=[]

        #Set correct values
        self.defineCorrectValues(exp_dict)

        self.myDetector.defineCorrectValuesDetector()
        self.mySource.defineCorrectValuesSource()
        self.mySampleofInterest.defineCorrectValuesSample()
        self.myAirVolume.defineCorrectValuesSample()
        self.myAirVolume.myThickness=(self.distSourceToMembrane+self.distObjectToDetector+self.distMembraneToObject)*1e6
        if self.myPlaque is not None:
            self.myPlaque.defineCorrectValuesSample()
        self.myMembrane.defineCorrectValuesSample()
        
        self.magnification=(self.distSourceToMembrane+self.distObjectToDetector+self.distMembraneToObject)/(self.distSourceToMembrane+self.distMembraneToObject)
        self.getStudyDimensions()
    
        #Set experiment data
        self.mySource.setMySpectrum()
        
        self.myAirVolume.getDeltaBeta(self.mySource.mySpectrum)
        self.myAirVolume.getMyGeometry(self.studyDimensions,self.studyPixelSize,self.sampling)
        if self.myPlaque is not None:
            self.myPlaque.getDeltaBeta(self.mySource.mySpectrum)
            self.myPlaque.getMyGeometry(self.studyDimensions,self.studyPixelSize,self.sampling)
        self.mySampleofInterest.getDeltaBeta(self.mySource.mySpectrum)
        self.mySampleofInterest.getMyGeometry(self.studyDimensions,self.studyPixelSize,self.sampling)

        self.myMembrane.getDeltaBeta(self.mySource.mySpectrum)
        self.myMembrane.membranePixelSize=self.studyPixelSize*self.distSourceToMembrane/(self.distSourceToMembrane+self.distMembraneToObject)
        self.myMembrane.getMyGeometry(self.studyDimensions,self.myMembrane.membranePixelSize,self.sampling, self.nbPoints)
        
        
        if self.nbPoints==0:
            print('Current experiment:',self.name)
            print("\nCurrent detector: ",self.myDetector.myName)
            print("  Detectors dimensions ",self.myDetector.myDimensions)
            print("Current source: ",self.mySource.myName)
            print("  Source type:",self.mySource.myType)
            print("  Source spectrum:",self.mySource.mySpectrum,"eV")
            print("Current sample:", self.mySampleofInterest.myName)
            print("  My sample type:", self.mySampleofInterest.myType)
            print("Current membrane:", self.myMembrane.myName)
            print("  My membrane type:", self.myMembrane.myType)
            print("Magnification :", self.magnification)
    
            print("Study dimensions:", self.studyDimensions)
            print("Study PixelSize =",self.studyPixelSize,"um")
            print("Over-sampling factor: ",self.sampling)
            
        
    def defineCorrectValues(self, exp_dict):
        #Initialize object source and detector
        self.mySource=Source()
        self.myDetector=Detector(exp_dict)
        
        for experiment in self.xmldoc.documentElement.getElementsByTagName("experiment"):
            correctExperiment = self.getText(experiment.getElementsByTagName("name")[0])
            if correctExperiment == self.name:
                self.distSourceToMembrane=float(self.getText(experiment.getElementsByTagName("distSourceToMembrane")[0]))
                self.distMembraneToObject=float(self.getText(experiment.getElementsByTagName("distMembraneToObject")[0]))
                self.distObjectToDetector=float(self.getText(experiment.getElementsByTagName("distObjectToDetector")[0]))
                self.meanShotCount=float(self.getText(experiment.getElementsByTagName("meanShotCount")[0]))/self.sampling**2
                
                for node in experiment.childNodes:
                    if node.localName=="plaqueName":
                        self.myPlaque=AnalyticalSample()
                        self.myPlaque.myName=self.getText(experiment.getElementsByTagName("plaqueName")[0])
                        
                self.myAirVolume=AnalyticalSample()
                self.myAirVolume.myName="air_volume"

                #Initializing object sample and membrane
                self.mySampleType=self.getText(experiment.getElementsByTagName("sampleType")[0])
                if self.mySampleType=="AnalyticalSample":
                    self.mySampleofInterest=AnalyticalSample()
                else:
                    raise Exception("sample type not defined")
                self.myMembrane=AnalyticalSample()
                
                #Getting the identity of the objects
                self.myMembrane.myName=self.getText(experiment.getElementsByTagName("membraneName")[0])
                self.mySampleofInterest.myName=self.getText(experiment.getElementsByTagName("sampleName")[0])
                self.myDetector.myName=self.getText(experiment.getElementsByTagName("detectorName")[0])
                self.mySource.myName=self.getText(experiment.getElementsByTagName("sourceName")[0])
                return
                        
        raise ValueError("experiment not found in xml file")
            
    
    def getText(self,node):
        return node.childNodes[0].nodeValue
    
    
    def getStudyDimensions(self):
        self.precision=(self.myDetector.myPixelSize/self.sampling/self.distObjectToDetector)
        self.studyDimensions=self.myDetector.myDimensions*int(self.sampling)
        self.studyDimensions[0]=int(self.studyDimensions[0])
        self.studyDimensions[1]=int(self.studyDimensions[1])
        self.studyPixelSize=self.myDetector.myPixelSize/self.sampling/self.magnification
        
    
    def wavePropagation(self, waveToPropagate, propagationDistance, Energy, magnification):
        if propagationDistance==0:
            return waveToPropagate
        
        #Propagateur de Fresnel
        Lambda=6.626*1e-34*2.998e8/(Energy*1000*1.6e-19)
        k=2*pi/Lambda
        
        Nx=self.studyDimensions[0]        
        Ny=self.studyDimensions[1]
        u, v = np.meshgrid(np.arange(0, Nx), np.arange(0, Ny))
        u = (u - (Nx / 2))
        v = (v - (Ny / 2))
        u_m = u *2*pi / (self.studyDimensions[0]*self.studyPixelSize*1e-6)
        v_m = v *2*pi / (self.studyDimensions[1]*self.studyPixelSize*1e-6)
        uv_sqr=  np.transpose(u_m ** 2 + v_m ** 2)  # ie (u2+v2)
        
        waveAfterPropagation=np.exp(1j*k*propagationDistance/magnification)*ifft2(ifftshift(np.exp(-1j*propagationDistance*(uv_sqr)/(2*k*magnification))*fftshift(fft2(waveToPropagate))))
        
        return waveAfterPropagation
    

    
    def refraction(self,intensityRefracted, phi, propagationDistance,Energy, magnification):
        timeBR=time.time()
        Lambda=6.626*1e-34*2.998e8/(Energy*1000*1.6e-19)
        k=2*pi/Lambda
        Nx,Ny=intensityRefracted.shape
        if propagationDistance==0:
            return intensityRefracted, 0, 0
        
        intensityRefracted2=np.zeros((intensityRefracted.shape))
        dphix,dphiy=np.gradient(phi,self.studyPixelSize*1e-6, edge_order=2)
        Dx=dphix*propagationDistance/k/(self.studyPixelSize*1e-6*magnification)
        Dy=dphiy*propagationDistance/k/(self.studyPixelSize*1e-6*magnification)        
        alphax=Dx
        alphay=Dy
        
        for i in range(Nx):
            for j in range(Ny):
                Dxtmp=Dx[i,j]
                Dytmp=Dy[i,j]
                inew=i
                jnew=j
                #Calculating displacement bigger than a pixel
                if abs(Dxtmp)>1:
                    inew=i+int(np.floor(Dxtmp))
                    Dxtmp=Dxtmp-np.floor(Dxtmp)
                if abs(Dytmp)>1:
                    jnew=j+int(np.floor(Dytmp))
                    Dytmp=Dytmp-np.floor(Dytmp)
                #Calculating sub-pixel displacement
                if 0<=inew<Nx:
                    if 0<=jnew<Ny:
                        intensityRefracted2[inew,jnew]+=intensityRefracted[i,j]*(1-abs(Dxtmp))*(1-abs(Dytmp))
                        
                        if inew<Nx-1:
                            if Dxtmp>=0:
                                if jnew<Ny-1:
                                    if Dytmp>=0:
                                        intensityRefracted2[inew+1,jnew]+=intensityRefracted[i,j]*Dxtmp*(1-Dytmp)
                                        intensityRefracted2[inew+1,jnew+1]+=intensityRefracted[i,j]*Dxtmp*Dytmp
                                        intensityRefracted2[inew,jnew+1]+=intensityRefracted[i,j]*(1-Dxtmp)*Dytmp
                        
                        if inew>0:
                            if Dxtmp<0:
                                if jnew<Ny-1:
                                    if Dytmp>=0:
                                        intensityRefracted2[inew-1,jnew]+=intensityRefracted[i,j]*abs(Dxtmp)*(1-Dytmp)
                                        intensityRefracted2[inew-1,jnew+1]+=intensityRefracted[i,j]*abs(Dxtmp)*Dytmp
                                        intensityRefracted2[inew,jnew+1]+=intensityRefracted[i,j]*(1-abs(Dxtmp))*Dytmp
                        
                        if inew>0:
                            if Dxtmp<0:
                                if jnew>0:
                                    if Dytmp<0:
                                        intensityRefracted2[inew-1,jnew]+=intensityRefracted[i,j]*abs(Dxtmp)*(1-abs(Dytmp))
                                        intensityRefracted2[inew-1,jnew-1]+=intensityRefracted[i,j]*Dxtmp*Dytmp
                                        intensityRefracted2[inew,jnew-1]+=intensityRefracted[i,j]*(1-abs(Dxtmp))*abs(Dytmp)
                        
                        if inew<Nx-1:
                            if Dxtmp>=0:
                                if jnew>0:
                                    if Dytmp<0:
                                        intensityRefracted2[inew+1,jnew]+=intensityRefracted[i,j]*Dxtmp*(1-abs(Dytmp))
                                        intensityRefracted2[inew+1,jnew-1]+=intensityRefracted[i,j]*Dxtmp*abs(Dytmp)
                                        intensityRefracted2[inew,jnew-1]+=intensityRefracted[i,j]*(1-Dxtmp)*abs(Dytmp)
        print("Refraction Time", time.time()-timeBR)
        return intensityRefracted2,alphax, alphay    
    
#        
    def computeSampleAndReferenceImages(self):
        
        #INITIALIZING PARAMETERS
        totalFlux=0 
        sumIntensity=0
        
        #INITIALIZING IMAGES
        incidentWave0=np.ones((self.studyDimensions[0],self.studyDimensions[1]))*np.sqrt(self.meanShotCount)
        self.imageSampleBeforeDetection=np.zeros((self.studyDimensions[0],self.studyDimensions[1]))
        self.imageReferenceBeforeDetection=np.zeros((self.studyDimensions[0],self.studyDimensions[1]))
        self.imagePropagBeforeDetection=np.zeros((self.studyDimensions[0],self.studyDimensions[1]))
        white=np.zeros((self.studyDimensions[0],self.studyDimensions[1]))
        
        SampleImage=np.zeros((self.myDetector.myDimensions[0]-2*self.myDetector.margins,self.myDetector.myDimensions[1]-2*self.myDetector.margins))
        ReferenceImage=np.zeros((self.myDetector.myDimensions[0]-2*self.myDetector.margins,self.myDetector.myDimensions[1]-2*self.myDetector.margins))
        PropagImage=np.zeros((self.myDetector.myDimensions[0]-2*self.myDetector.margins,self.myDetector.myDimensions[1]-2*self.myDetector.margins))
        detectedWhite=np.zeros((self.myDetector.myDimensions[0]-2*self.myDetector.margins,self.myDetector.myDimensions[1]-2*self.myDetector.margins))
        
        
        #Defining total flux for normalizing spectrum
        for currentEnergy, flux in self.mySource.mySpectrum:
            totalFlux+=flux
        i=0
        #Calculating everything for each energy of the spectrum
        for currentEnergy, flux in self.mySource.mySpectrum:
            if currentEnergy<=self.myDetector.myEnergyLimit:
                #Taking into account source window and air attenuation of intensity
                incidentIntensity=incidentWave0**2
                incidentIntensity=incidentIntensity*flux/totalFlux
                incidentIntensity, _=self.myAirVolume.setWaveRT(incidentIntensity,1, currentEnergy)
                incidentWave=np.sqrt(incidentIntensity)
                
                print("\nCurrent Energy:", currentEnergy)
                
                #Passage of the incident wave through the membrane
                print("Setting wave through membrane")
                self.waveSampleAfterMembrane=self.myMembrane.setWave(incidentWave,currentEnergy)
                            
                magMemObj=(self.distSourceToMembrane+self.distMembraneToObject)/self.distSourceToMembrane
                self.waveSampleBeforeSample=self.wavePropagation(self.waveSampleAfterMembrane,self.distMembraneToObject,currentEnergy,magMemObj)
    
                print("Setting wave through sample for sample image")            
                self.waveSampleAfterSample=self.mySampleofInterest.setWave(self.waveSampleBeforeSample,currentEnergy)
                
                #Propagation to detector
                print("Propagating waves to detector plane")
                self.waveSampleBeforeDetection=self.wavePropagation(self.waveSampleAfterSample,self.distObjectToDetector,currentEnergy,self.magnification)
                self.waveReferenceBeforeDetection=self.wavePropagation(self.waveSampleAfterMembrane,self.distObjectToDetector+self.distMembraneToObject,currentEnergy,self.magnification)
                #Combining intensities for several energies
                intensitySampleBeforeDetection=abs(self.waveSampleBeforeDetection)**2
                if self.myPlaque is not None:
                    intensitySampleBeforeDetection,_=self.myPlaque.setWaveRT(intensitySampleBeforeDetection,1, currentEnergy)
                intensityReferenceBeforeDetection=abs(self.waveReferenceBeforeDetection)**2
                if self.myPlaque is not None:
                    intensityReferenceBeforeDetection,_=self.myPlaque.setWaveRT(intensityReferenceBeforeDetection,1, currentEnergy)
                self.imageSampleBeforeDetection+=intensitySampleBeforeDetection
                self.imageReferenceBeforeDetection+=intensityReferenceBeforeDetection
                
                sumIntensity+=np.mean(intensityReferenceBeforeDetection)
                self.meanEnergy+=currentEnergy*np.mean(intensityReferenceBeforeDetection)
                
                if self.nbPoints==0: #We only do it for the first point
                    print("Setting wave through sample for propag and abs image")
                    self.wavePropagAfterSample=self.mySampleofInterest.setWave(incidentWave,currentEnergy)
                    self.wavePropagBeforeDetection=self.wavePropagation(self.wavePropagAfterSample,self.distObjectToDetector,currentEnergy,self.magnification)
                    intensityPropagBeforeDetection=abs(self.wavePropagBeforeDetection)**2
                    if self.myPlaque is not None:
                        intensityPropagBeforeDetection,_=self.myPlaque.setWaveRT(intensityPropagBeforeDetection,1, currentEnergy)
                    self.imagePropagBeforeDetection+=intensityPropagBeforeDetection
                    i+=1
                    incidentIntensityWhite=incidentWave**2
                    if self.myPlaque is not None:
                        incidentIntensityWhite,_=self.myPlaque.setWaveRT(incidentWave**2,1, currentEnergy)
                    white+=incidentIntensityWhite
                
                    
        self.meanEnergy=self.meanEnergy/sumIntensity
        print("Mean energy detected in reference image", self.meanEnergy)
            
        effectiveSourceSize=self.mySource.mySize*self.distObjectToDetector/(self.distSourceToMembrane+self.distMembraneToObject)/self.myDetector.myPixelSize*self.sampling #FWHM
        self.imageSampleBeforeDetection=self.imageSampleBeforeDetection
        self.imageReferenceBeforeDetection=self.imageReferenceBeforeDetection
        imageSampleBeforeDetection=self.imageSampleBeforeDetection
        #DETECTION IMAGES FOR ENERGY BIN
        print("Detection sample image")
        SampleImage=self.myDetector.detection(self.imageSampleBeforeDetection,effectiveSourceSize)
        print("Detection reference image")
        ReferenceImage=self.myDetector.detection(self.imageReferenceBeforeDetection,effectiveSourceSize)
        if self.nbPoints==0:
            self.imagePropagBeforeDetection=self.imagePropagBeforeDetection
            print("Detection propagation image")
            PropagImage=self.myDetector.detection(self.imagePropagBeforeDetection,effectiveSourceSize)
            detectedWhite=self.myDetector.detection(white,effectiveSourceSize)

        return  SampleImage, ReferenceImage,PropagImage,detectedWhite
        
    def computeSampleAndReferenceImagesRT(self):
        
        #INITIALIZING PARAMETERS
        totalFlux=0 
        sumIntensity=0
        
        #INITIALIZING IMAGES
        incidentIntensity0=np.ones((self.studyDimensions[0],self.studyDimensions[1]))*(self.meanShotCount)
        incidentPhi=np.zeros((self.studyDimensions[0],self.studyDimensions[1]))
        self.imageSampleBeforeDetection=np.zeros((self.studyDimensions[0],self.studyDimensions[1]))
        self.imageReferenceBeforeDetection=np.zeros((self.studyDimensions[0],self.studyDimensions[1]))
        self.imagePropagBeforeDetection=np.zeros((self.studyDimensions[0],self.studyDimensions[1]))
        white=np.zeros((self.studyDimensions[0],self.studyDimensions[1]))
        
        SampleImage=np.zeros((self.myDetector.myDimensions[0]-2*self.myDetector.margins,self.myDetector.myDimensions[1]-2*self.myDetector.margins))
        ReferenceImage=np.zeros((self.myDetector.myDimensions[0]-2*self.myDetector.margins,self.myDetector.myDimensions[1]-2*self.myDetector.margins))
        PropagImage=np.zeros((self.myDetector.myDimensions[0]-2*self.myDetector.margins,self.myDetector.myDimensions[1]-2*self.myDetector.margins))
        detectedWhite=np.zeros((self.myDetector.myDimensions[0]-2*self.myDetector.margins,self.myDetector.myDimensions[1]-2*self.myDetector.margins))

        #Defining total flux for normalizing spectrum
        for currentEnergy, flux in self.mySource.mySpectrum:
            totalFlux+=flux
            
        #Calculating everything for each energy of the spectrum
        for currentEnergy, flux in self.mySource.mySpectrum:
            if currentEnergy<=self.myDetector.myEnergyLimit:
                incidentIntensity=incidentIntensity0*flux/totalFlux
                incidentIntensity, _=self.myAirVolume.setWaveRT(incidentIntensity,1, currentEnergy)
                
                print("\nCurrent Energy:", currentEnergy)
                
                #Passage of the incident wave through the membrane
                print("Setting wave through membrane")
                self.IntensitySampleAfterMembrane, phiWaveSampleAfterMembrane=self.myMembrane.setWaveRT(incidentIntensity,incidentPhi,currentEnergy)
                
                #Propagation from membrane to object and passage through the object
                self.IntensitySampleBeforeSample,_,_=self.refraction(abs(self.IntensitySampleAfterMembrane),phiWaveSampleAfterMembrane,self.distMembraneToObject,currentEnergy,self.magnification)
    
                print("Setting wave through sample for sample image")            
                self.IntensitySampleAfterSample,phiWaveSampleAfterSample=self.mySampleofInterest.setWaveRT(self.IntensitySampleBeforeSample,phiWaveSampleAfterMembrane,currentEnergy)
                
                #Propagation to detector
                print("Propagating waves to detector plane")
                self.imageSampleAfterRefraction,_,_=self.refraction(abs(self.IntensitySampleAfterSample),phiWaveSampleAfterSample,self.distObjectToDetector,currentEnergy,self.magnification)
                self.imageReferenceAfterRefraction,_,_=self.refraction(abs(self.IntensitySampleBeforeSample),phiWaveSampleAfterMembrane,self.distObjectToDetector,currentEnergy,self.magnification)
                intensitySampleBeforeDetection=self.imageSampleAfterRefraction
                intensityReferenceBeforeDetection=self.imageReferenceAfterRefraction
                #Plaque attenuation
                if self.myPlaque is not None:
                    intensitySampleBeforeDetection,_=self.myPlaque.setWaveRT(intensitySampleBeforeDetection,1, currentEnergy)
                    intensityReferenceBeforeDetection,_=self.myPlaque.setWaveRT(intensityReferenceBeforeDetection,1, currentEnergy)
                #Combining intensities for several energies
                self.imageSampleBeforeDetection+=intensitySampleBeforeDetection
                self.imageReferenceBeforeDetection+=intensityReferenceBeforeDetection
                
                sumIntensity+=np.mean(intensityReferenceBeforeDetection)
                self.meanEnergy+=currentEnergy*np.mean(intensityReferenceBeforeDetection)
                
                if self.nbPoints==0: #We only do it for the first point
                    print("Setting wave through sample for propag and abs image")
                    self.IntensityPropagAfterSample,phiWavePropagAfterSample=self.mySampleofInterest.setWaveRT(incidentIntensity,incidentPhi,currentEnergy)
                    self.imagePropagAfterRefraction, self.Dxreal, self.Dyreal=self.refraction(abs(self.IntensityPropagAfterSample),phiWavePropagAfterSample,self.distObjectToDetector,currentEnergy,self.magnification)
                    intensityPropagBeforeDetection=self.imagePropagAfterRefraction
                    if self.myPlaque is not None:
                        intensityPropagBeforeDetection,_=self.myPlaque.setWaveRT(intensityPropagBeforeDetection,1, currentEnergy)
                        incidentIntensity,_=self.myPlaque.setWaveRT(incidentIntensity,1, currentEnergy)
                    white+=incidentIntensity
                    self.imagePropagBeforeDetection+=intensityPropagBeforeDetection



        self.meanEnergy=self.meanEnergy/sumIntensity
                
        self.meanEnergy=self.meanEnergy/totalFlux
        print("MeanEnergy", self.meanEnergy)
        
        effectiveSourceSize=self.mySource.mySize*self.distObjectToDetector/(self.distSourceToMembrane+self.distMembraneToObject)/self.myDetector.myPixelSize*self.sampling #FWHM
        self.imageSampleBeforeDetection=self.imageSampleBeforeDetection
        self.imageReferenceBeforeDetection=self.imageReferenceBeforeDetection

        #DETECTION IMAGES FOR ENERGY BIN
        print("Detection sample image")
        SampleImage=self.myDetector.detection(self.imageSampleBeforeDetection,effectiveSourceSize)
        print("Detection reference image")
        ReferenceImage=self.myDetector.detection(self.imageReferenceBeforeDetection,effectiveSourceSize)
        
        if self.nbPoints==0:
            self.imagePropagBeforeDetection=self.imagePropagBeforeDetection
            print("Detection propagation image")
            PropagImage=self.myDetector.detection(self.imagePropagBeforeDetection,effectiveSourceSize)
                
                
        detectedWhite=self.myDetector.detection(white,effectiveSourceSize)
                
        return  SampleImage, ReferenceImage,PropagImage,detectedWhite, self.Dxreal, self.Dyreal

    
        
    def saveAllParameters(self,time0,expDict):
        fileName=expDict['filepath']+self.name+'_'+str(expDict['expID'])+".txt"
        print("file name: ", fileName)
        f=open(fileName,"w+")

        f.write("EXPERIMENT PARAMETERS - Fresnel - "+str(expDict['expID']))

        f.write("\n\nDistances:")
        f.write("\nDistance source to membrane: %gm" %self.distSourceToMembrane)
        f.write("\nDistance membrane to sample: %gm" %self.distMembraneToObject)
        f.write("\nDistance sample to detector: %gm" %self.distObjectToDetector)

        f.write("\n\nSource parameters:")
        f.write("\nSource name: %s" %self.mySource.myName)
        f.write("\nSource type: %s" %self.mySource.myType)
        f.write("\nSource size: %gum" %self.mySource.mySize)
        if self.mySource.myType=="Monochromatic":
            f.write("\nSource energy: %gkev" %(self.mySource.mySpectrum[0][0]/1000))
        if self.mySource.myType=="Polychromatic":
            f.write("\nSource voltage: %gkVp" %self.mySource.myVoltage)
            f.write("\nSource spectrum energy sampling: %gkeV" %self.mySource.myEnergySampling)
        
        f.write("\n\nDetector parameters:")
        f.write("\nDetector name: %s"%self.myDetector.myName)
        f.write("\nDetector dimensions:"+str(self.myDetector.myDimensions[0]-expDict['margin'])+"x"+str(self.myDetector.myDimensions[1]-expDict['margin'])+"pix")
        f.write("\nDetector pixels size: %gum" %self.myDetector.myPixelSize)
        f.write("\nDetector PSF: %gpix" %self.myDetector.myPSF)
        
        f.write("\n\nSample informations")
        f.write("\nSample name: %s" %self.mySampleofInterest.myName)
        f.write("\nSample type: %s" %self.mySampleType)
        if self.mySampleofInterest.myGeometryFunction=="CreateSampleCylindre":
            f.write("\nWire's radius: %s" %self.mySampleofInterest.myRadius)
            f.write("\nWire's material: %s" %self.mySampleofInterest.myMaterials)
        if self.mySampleofInterest.myGeometryFunction=="openContrastPhantom":
            f.write("\nContrast Phantom geometry folder: %s" %self.mySampleofInterest.myGeometryFolder)
        
        f.write("\n\nMembrane informations:")
        f.write("\nMembrane name: %s" %self.myMembrane.myName)
        f.write("\nMembrane type: %s" %self.myMembrane.myType)
        f.write("\nMembrane geometry function: %s" %self.myMembrane.myGeometryFunction)
        if self.myMembrane.myGeometryFunction=="getMembraneFromFile":
            f.write("\nMembrane geometry file: %s" %self.myMembrane.myMembraneFile)
        if self.myMembrane.myGeometryFunction=="getMembraneSegmentedFromFile":
            f.write("\nMembrane number of layers: %s" %self.myMembrane.myNbOfLayers)  
            
        if self.myPlaque is not None:
            f.write("\n\nDetectors protection plaque")
            f.write("Plaque thickness: %s" %self.myPlaque.myThickness)
            f.write("Plaque Material: %s" %self.myPlaque.myMaterials)

            
        f.write("\n\nOther study parameters:")
        f.write("\nPrecision of the calculation: %gurad" %self.precision)
        f.write("\nOversampling Factor: %g" %self.sampling)
        f.write("\nStudy dimensions: "+str(self.studyDimensions[0])+"x"+str(self.studyDimensions[1])+"pix")
        f.write("\nStudy pixel size: %gum" %self.studyPixelSize)
        f.write("\nMean shot count: %g" %self.meanShotCount*self.sampling**2)        
        f.write("\nNumber of points: %g" %(self.nbPoints+1))       
        f.write("\nEntire computing time: %gs" %(time.time()-time0))   
        
        f.close()
        
        return
    
    
    def createExpDict(self, expParam):
        expParam['energy'] = self.meanEnergy
        expParam['pixel'] = self.myDetector.myPixelSize * 1e-6
        expParam['distOD'] = self.distObjectToDetector
        expParam['distSO'] = self.distSourceToMembrane+self.distMembraneToObject
        expParam['sourceSize'] = self.mySource.mySize * 1e-6  # in meters
        expParam['detectorPSF'] = self.myDetector.myPSF  # in pixel
    
        return expParam
    
        