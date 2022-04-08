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
from usefullScripts.getSamplingFactor import is_overSampling_ok
from refractionFileNumba2 import fastRefraction, fastRefractionDF
from getk import getk

class Experiment:
    def __init__(self, exp_dict):#expName, pointNum, overSampling):
        """Experiment class constructor.
            
        Args:
            exp_dict (dict): dictionnary of simulation algorithm parameters
        """
        self.xmlExperimentFileName="xmlFiles/Experiment.xml"
        self.xmldoc = minidom.parse(self.xmlExperimentFileName)
        self.name=exp_dict['experimentName']
        self.exp_dict=exp_dict
        self.exp_dict['studyPixelSize']=0.
        self.exp_dict['studyDimensions']=(0.,0.)
        self.exp_dict['inVacuum']=False
        self.exp_dict['meanShotCount']=0
        self.exp_dict['meanEnergy']=0
        self.exp_dict['distSourceToMembrane']=0
        self.exp_dict['distMembraneToObject']=0
        self.exp_dict['distObjectToDetector']=0
        
        #Parameters units
        self.exp_dict['studyPixelSize_unit']="um"
        self.exp_dict['studyDimensions_unit']="pixels"
        self.exp_dict['meanEnergy_unit']="keV"
        self.exp_dict['distSourceToMembrane_unit']="m"
        self.exp_dict['distMembraneToObject_unit']="m"
        self.exp_dict['distObjectToDetector_unit']="m"
        
        
        self.mySampleofInterest=None
        self.mySampleType=""
        self.myDetector=None
        self.mySource=None
        self.myMembrane=None
        self.myPlate=None
        self.myAirVolume=None
        
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
        self.myAirVolume.myThickness=(self.exp_dict['distSourceToMembrane']+self.exp_dict['distObjectToDetector']+self.exp_dict['distMembraneToObject'])*1e6
        if self.myPlate is not None:
            self.myPlate.defineCorrectValuesSample()
        self.myMembrane.defineCorrectValuesSample()
        
        self.exp_dict['magnification']=(self.exp_dict['distSourceToMembrane']+self.exp_dict['distObjectToDetector']+self.exp_dict['distMembraneToObject'])/(self.exp_dict['distSourceToMembrane']+self.exp_dict['distMembraneToObject'])
        self.getStudyDimensions()
    
        #Set experiment data
        self.mySource.setMySpectrum(self.myDetector.det_param["photonCounting"])
        
        self.myAirVolume.getDeltaBeta(self.mySource.mySpectrum)
        self.myAirVolume.getMyGeometry(self.exp_dict['studyDimensions'],self.exp_dict['studyPixelSize'],self.exp_dict['overSampling'])
        if self.myPlate is not None:
            self.myPlate.getDeltaBeta(self.mySource.mySpectrum)
            self.myPlate.getMyGeometry(self.exp_dict['studyDimensions'],self.exp_dict['studyPixelSize'],self.exp_dict['overSampling'])
        self.mySampleofInterest.getDeltaBeta(self.mySource.mySpectrum)
        self.mySampleofInterest.getMyGeometry(self.exp_dict['studyDimensions'],self.exp_dict['studyPixelSize'],self.exp_dict['overSampling'])

        self.myMembrane.getDeltaBeta(self.mySource.mySpectrum)
        self.myMembrane.membranePixelSize=self.exp_dict['studyPixelSize']*self.exp_dict['distSourceToMembrane']/(self.exp_dict['distSourceToMembrane']+self.exp_dict['distMembraneToObject'])
        
        if self.myDetector.det_param['myScintillatorMaterial'] is not None:
            self.myDetector.getBeta(self.mySource.mySpectrum)
            self.myDetector.getSpectralEfficiency()
        
        #Check if oversampling is ok
        if self.exp_dict['simulation_type']=="RayT":
            if self.exp_dict["overSampling"]<2:
                print(f'/!\/!\ OVERSAMPLING FACTOR < MIN OVERSAMPLING FOR RAY-T MODEL: {exp_dict["overSampling"]} < 2')
        if self.exp_dict['simulation_type']=="Fresnel":
            if self.mySource.source_dict["myType"]=="Polychromatic":
                is_overSampling_ok(self.exp_dict, self.myDetector.det_param['myPixelSize'], self.mySource.mySpectrum[-1][0]/2)
            elif self.mySource.source_dict["myType"]=="Monochromatic":
                is_overSampling_ok(self.exp_dict, self.myDetector.det_param['myPixelSize'], self.mySource.source_dict["Energy"])                
        
        
        print('\nCurrent experiment:',self.name)
        print(f'  Experiment in Vacuum: {self.exp_dict["inVacuum"]}')
        print("  Magnification :", self.exp_dict['magnification'])
        print(f'  Study dimensions: {self.exp_dict["studyDimensions"]} pixels')
        print("  Sample pixel size =",self.exp_dict["studyPixelSize"],"um")
        print("  Over-overSampling factor: ",self.exp_dict["overSampling"])
        
        print("\nCurrent detector: ",self.myDetector.myName)
        print(f'  Detector pixel size: {self.myDetector.det_param["myPixelSize"]} um')
        if self.myDetector.det_param['myScintillatorMaterial'] is not None:
            print(f'  Scintillator {self.myDetector.det_param["myScintillatorMaterial"]} of {self.myDetector.det_param["myScintillatorThickness"]}um')
        print("  Detectors dimensions: ",self.myDetector.det_param["myDimensions"])
        print("\nCurrent source: ",self.mySource.myName)
        print("  Source type:",self.mySource.source_dict["myType"])
        if self.mySource.source_dict["myType"]=='Monochromatic':
            print(f'Energy: {self.mySource.source_dict["Energy"]} keV')
        else:
            if "myVoltage" in self.mySource.source_dict:
                print(f'Source voltage: {self.mySource.source_dict["myVoltage"]} kVp')
            if "myTargetMaterial" in self.mySource.source_dict:
                print(f'Anode material: {self.mySource.source_dict["myTargetMaterial"]}')
            if self.mySource.source_dict["filterMaterial"] is not None:
                print(f'filter: {self.mySource.source_dict["filterMaterial"]} of {self.mySource.source_dict["filterThickness"]} mm')
        print("\nCurrent sample:", self.mySampleofInterest.myName)
        print("\nCurrent membrane:", self.myMembrane.myName)
            
        
    def defineCorrectValues(self, exp_dict):
        """
        Initializes every compound parameters before calculations

        Args:
            exp_dict (dictionnary): algorithm parameters.

        Raises:
            Exception: sample type not defined.
            ValueError: experiment not found in xml file.

        Returns:
            None.

        """
        #Initialize object source and detector
        self.mySource=Source()
        self.myDetector=Detector(exp_dict)
        
        for experiment in self.xmldoc.documentElement.getElementsByTagName("experiment"):
            correctExperiment = self.getText(experiment.getElementsByTagName("name")[0])
            if correctExperiment == self.name:
                self.exp_dict['distSourceToMembrane']=float(self.getText(experiment.getElementsByTagName("distSourceToMembrane")[0]))
                self.exp_dict['distMembraneToObject']=float(self.getText(experiment.getElementsByTagName("distMembraneToObject")[0]))
                self.exp_dict['distObjectToDetector']=float(self.getText(experiment.getElementsByTagName("distObjectToDetector")[0]))
                self.exp_dict['meanShotCount']=float(self.getText(experiment.getElementsByTagName("meanShotCount")[0]))#/self.exp_dict['overSampling']**2
                
                for node in experiment.childNodes:
                    if node.localName=="inVacuum":
                        text=self.getText(experiment.getElementsByTagName("inVacuum")[0])
                        if text=="True":
                            self.exp_dict['inVacuum']=True
                        else:
                            self.exp_dict['inVacuum']=False
                    if node.localName=="plateName":
                        self.myPlate=AnalyticalSample()
                        self.myPlate.myName=self.getText(experiment.getElementsByTagName("plateName")[0])
                        
                        
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
        """
        Calculates the study dimensions considereing the geometry of the set up, the field of view and the sample pixels overoverSampling

        Returns:
            None.

        """
        self.precision=(self.myDetector.det_param["myPixelSize"]/self.exp_dict['overSampling']/self.exp_dict['distObjectToDetector'])
        self.exp_dict['studyDimensions']=self.myDetector.det_param["myDimensions"]*int(self.exp_dict['overSampling'])
        self.exp_dict['studyDimensions'][0]=int(self.exp_dict['studyDimensions'][0])
        self.exp_dict['studyDimensions'][1]=int(self.exp_dict['studyDimensions'][1])
        self.exp_dict['studyPixelSize']=self.myDetector.det_param["myPixelSize"]/self.exp_dict['overSampling']/self.exp_dict['magnification']
        
    
    def wavePropagation(self, waveToPropagate, propagationDistance, Energy, magnification):
        """
        Propagation of the wave 

        Args:
            waveToPropagate (2d numpy array): incident wave.
            propagationDistance (Float): propagation distance in m.
            Energy (Float): considered eneregy in keV.
            magnification (Float): magnification on the considered segment from the source.

        Returns:
            TYPE: DESCRIPTION.

        """
        if propagationDistance==0:
            return waveToPropagate
        
        margin=15
        waveToPropagate=np.pad(waveToPropagate,margin, mode='reflect')
        #Propagateur de Fresnel
        k=getk(Energy*1000)
        Nx,Ny=waveToPropagate.shape
        # Nx=self.exp_dict['studyDimensions'][0]        
        # Ny=self.exp_dict['studyDimensions'][1]
        u, v = np.meshgrid(np.arange(0, Nx), np.arange(0, Ny))
        u = (u - (Nx / 2))
        v = (v - (Ny / 2))
        u_m = u *2*pi / (self.exp_dict['studyDimensions'][0]*self.exp_dict['studyPixelSize']*1e-6)
        v_m = v *2*pi / (self.exp_dict['studyDimensions'][1]*self.exp_dict['studyPixelSize']*1e-6)
        uv_sqr=  np.transpose(u_m ** 2 + v_m ** 2)  # ie (u2+v2)
        
        waveAfterPropagation=np.exp(1j*k*propagationDistance/magnification)*ifft2(ifftshift(np.exp(-1j*propagationDistance*(uv_sqr)/(2*k*magnification))*fftshift(fft2(waveToPropagate))))
        waveAfterPropagation=waveAfterPropagation[margin:Nx-margin,margin:Ny-margin]
        return waveAfterPropagation
    
    
    def refraction(self,intensityRefracted, phi, propagationDistance,Energy, magnification, darkField=0):
        """
        Computes the intensity after propagation with ray-tracing

        Args:
            intensityRefracted (2d numpy array): intensity before propagation.
            phi (2d numpy array): phase.
            propagationDistance (Float): propagation distance in m.
            Energy (Float): considered energy of the spectrum (in keV).
            magnification (Float): magnification of the considered segment from the source.

        Returns:
            intensityRefracted2 (2d numpy array): DESCRIPTION.
            Dx (2d numpy array): displacements along x.
            Dy (2d numpy array): displacements along y.

        """
        if type(darkField) == int or type(darkField)==float:
            intensityRefracted2,Dx, Dy=fastRefraction(intensityRefracted, phi, propagationDistance,Energy, magnification,self.exp_dict["studyPixelSize"])
        else:
            intensityRefracted2,Dx, Dy=fastRefractionDF(intensityRefracted, phi, propagationDistance,Energy, magnification,self.exp_dict["studyPixelSize"], darkField)
        
        return intensityRefracted2,Dx, Dy   
#        
    def computeSampleAndReferenceImages_Fresnel(self, pointNum):
        """
        Compute intensity changes on the path of the previously difined experiment 
        to create all the images of the SBI experiment with FRESNEL PROPAGATOR


        Returns:
            SampleImage (2d numpy array): sample image simulated with sample and membrane.
            ReferenceImage (2d numpy array): reference image simulated with membrane.
            PropagImage (2d numpy array): propagation image simulated with only sample.
            detectedWhite (2d numpy array): white image without membrane and sample.

        """
        
        #INITIALIZING PARAMETERS
        sumIntensity=0
        ibin=0
        if pointNum==0:
            if any(elem<self.mySource.mySpectrum[0][0] for elem in self.myDetector.det_param["myBinsThersholds"]) or any(elem>self.mySource.mySpectrum[-1][0] for elem in self.myDetector.det_param["myBinsThersholds"]):
                raise Exception(f'At least one of your detector bin threshold is outside your source spectrum. \nYour source spectrum ranges from {self.mySource.mySpectrum[0][0]} to {self.mySource.mySpectrum[-1][0]}')
            # self.myDetector.det_param["myBinsThersholds"].insert(0,self.mySource.mySpectrum[0][0])
            self.myDetector.det_param["myBinsThersholds"].append(self.mySource.mySpectrum[-1][0])
        nbins=len(self.myDetector.det_param["myBinsThersholds"])
        SampleImage=np.zeros((nbins,self.myDetector.det_param['myDimensions'][0],self.myDetector.det_param['myDimensions'][1]))
        ReferenceImage=np.zeros((nbins,self.myDetector.det_param['myDimensions'][0],self.myDetector.det_param['myDimensions'][1]))
        PropagImage=np.zeros((nbins, self.myDetector.det_param['myDimensions'][0],self.myDetector.det_param['myDimensions'][1]))
        detectedWhite=np.zeros((nbins, self.myDetector.det_param['myDimensions'][0],self.myDetector.det_param['myDimensions'][1]))
        
        #INITIALIZING IMAGES
        incidentIntensity0=np.ones((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))*(self.exp_dict['meanShotCount']/self.exp_dict['overSampling']**2)
        self.imageSampleBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        self.imageReferenceBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        self.imagePropagBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        white=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        
        
        i=0
        #Calculating everything for each energy of the spectrum
        for currentEnergy, flux in self.mySource.mySpectrum:
            print("Current Energy:", currentEnergy)
            #Taking into account source window and air attenuation of intensity
            incidentIntensity=incidentIntensity0*flux
            
            if not self.exp_dict['inVacuum']:
                incidentIntensity, _,_=self.myAirVolume.setWaveRT(incidentIntensity,currentEnergy)
            
            #Take into account the detector scintillator efficiency if given in xml file
            if self.myDetector.det_param['myScintillatorMaterial'] is not None:
                for energyData, betaEn in self.myDetector.beta:
                    if energyData==currentEnergy:
                        beta=betaEn
                k=getk(currentEnergy*1000)
                detectedSpectrum=1-np.exp(-2*k*self.myDetector.det_param['myScintillatorThickness']*1e-6*beta)
                # print("Scintillator efficiency for current energy:", detectedSpectrum)
                incidentIntensity=incidentIntensity*detectedSpectrum
            incidentWave=np.sqrt(incidentIntensity)
            
            #Passage of the incident wave through the membrane
            # print("Setting wave through membrane")
            self.waveSampleAfterMembrane=self.myMembrane.setWave(incidentWave,currentEnergy)
                        
            magMemObj=(self.exp_dict['distSourceToMembrane']+self.exp_dict['distMembraneToObject'])/self.exp_dict['distSourceToMembrane']
            self.waveSampleBeforeSample=self.wavePropagation(self.waveSampleAfterMembrane,self.exp_dict['distMembraneToObject'],currentEnergy,magMemObj)

            # print("Setting wave through sample for sample image")            
            self.waveSampleAfterSample=self.mySampleofInterest.setWave(self.waveSampleBeforeSample,currentEnergy)
            
            #Propagation to detector
            # print("Propagating waves to detector plane")
            self.waveSampleBeforeDetection=self.wavePropagation(self.waveSampleAfterSample,self.exp_dict['distObjectToDetector'],currentEnergy,self.exp_dict['magnification'])
            self.waveReferenceBeforeDetection=self.wavePropagation(self.waveSampleAfterMembrane,self.exp_dict['distObjectToDetector']+self.exp_dict['distMembraneToObject'],currentEnergy,self.exp_dict['magnification'])
            #Combining intensities for several energies
            intensitySampleBeforeDetection=abs(self.waveSampleBeforeDetection)**2
            if self.myPlate is not None:
                intensitySampleBeforeDetection,_,_=self.myPlate.setWaveRT(intensitySampleBeforeDetection,currentEnergy)
            intensityReferenceBeforeDetection=abs(self.waveReferenceBeforeDetection)**2
            if self.myPlate is not None:
                intensityReferenceBeforeDetection,_,_=self.myPlate.setWaveRT(intensityReferenceBeforeDetection, currentEnergy)
            self.imageSampleBeforeDetection+=intensitySampleBeforeDetection
            self.imageReferenceBeforeDetection+=intensityReferenceBeforeDetection
            
            sumIntensity+=np.mean(intensityReferenceBeforeDetection)
            self.exp_dict['meanEnergy']+=currentEnergy*np.mean(intensityReferenceBeforeDetection)
            
            if pointNum==0: #We only do it for the first point
                # print("Setting wave through sample for propag and abs image")
                self.wavePropagAfterSample=self.mySampleofInterest.setWave(incidentWave,currentEnergy)
                self.wavePropagBeforeDetection=self.wavePropagation(self.wavePropagAfterSample,self.exp_dict['distObjectToDetector'],currentEnergy,self.exp_dict['magnification'])
                intensityPropagBeforeDetection=abs(self.wavePropagBeforeDetection)**2
                if self.myPlate is not None:
                    intensityPropagBeforeDetection,_, _=self.myPlate.setWaveRT(intensityPropagBeforeDetection,currentEnergy)
                self.imagePropagBeforeDetection+=intensityPropagBeforeDetection
                i+=1
                incidentIntensityWhite=incidentWave**2
                if self.myPlate is not None:
                    incidentIntensityWhite,_, _=self.myPlate.setWaveRT(incidentWave**2, currentEnergy)
                white+=incidentIntensityWhite
                
                
            if currentEnergy>self.myDetector.det_param["myBinsThersholds"][ibin]-self.mySource.source_dict["myEnergySampling"]/2:
            
                effectiveSourceSize=self.mySource.source_dict["mySize"]*self.exp_dict['distObjectToDetector']/(self.exp_dict['distSourceToMembrane']+self.exp_dict['distMembraneToObject'])/self.myDetector.det_param['myPixelSize']*self.exp_dict['overSampling'] #FWHM
        
                #DETECTION IMAGES FOR ENERGY BIN
                    
                print("Mean energy detected in reference image", self.exp_dict['meanEnergy'])
                    
                #DETECTION IMAGES FOR ENERGY BIN
                # print("Detection sample image")
                SampleImage[ibin]=self.myDetector.detection(self.imageSampleBeforeDetection,effectiveSourceSize,self.exp_dict)
                # print("Detection reference image")
                ReferenceImage[ibin]=self.myDetector.detection(self.imageReferenceBeforeDetection,effectiveSourceSize,self.exp_dict)
                if pointNum==0:
                    # print("Detection propagation image")
                    PropagImage[ibin]=self.myDetector.detection(self.imagePropagBeforeDetection,effectiveSourceSize,self.exp_dict)
                detectedWhite[ibin]=self.myDetector.detection(white,effectiveSourceSize, self.exp_dict)
                
                self.imageSampleBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
                self.imageReferenceBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
                self.imagePropagBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
                white=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
                
                ibin+=1
            
        self.exp_dict['meanEnergy']=self.exp_dict['meanEnergy']/sumIntensity

        return  SampleImage, ReferenceImage,PropagImage,detectedWhite
        
    def computeSampleAndReferenceImages_RT(self, pointNum):
        """
        Compute intensity changes on the path of the previously difined experiment 
        to create all the images of the SBI experiment with ray-tracing

        Returns:
            SampleImage (2d numpy array): sample image simulated with sample and membrane.
            ReferenceImage (2d numpy array): reference image simulated with membrane.
            PropagImage (2d numpy array): propagation image simulated with only sample.
            detectedWhite (2d numpy array): white image without membrane and sample.
            2d numpy array: real Dx from sample to detector.
            2d numpy array: real Dy from sample to detector.

        """
        
        #INITIALIZING PARAMETERS
        sumIntensity=0
        ibin=0
        if pointNum==0:
            if any(elem<self.mySource.mySpectrum[0][0] for elem in self.myDetector.det_param["myBinsThersholds"]) or any(elem>self.mySource.mySpectrum[-1][0] for elem in self.myDetector.det_param["myBinsThersholds"]):
                raise Exception(f'At least one of your detector bin threshold is outside your source spectrum. \nYour source spectrum ranges from {self.mySource.mySpectrum[0][0]} to {self.mySource.mySpectrum[-1][0]}')
            # self.myDetector.det_param["myBinsThersholds"].insert(0,self.mySource.mySpectrum[0][0])
            self.myDetector.det_param["myBinsThersholds"].append(self.mySource.mySpectrum[-1][0])
        nbins=len(self.myDetector.det_param["myBinsThersholds"])
        SampleImage=np.zeros((nbins,self.myDetector.det_param['myDimensions'][0],self.myDetector.det_param['myDimensions'][1]))
        ReferenceImage=np.zeros((nbins,self.myDetector.det_param['myDimensions'][0],self.myDetector.det_param['myDimensions'][1]))
        PropagImage=np.zeros((nbins, self.myDetector.det_param['myDimensions'][0],self.myDetector.det_param['myDimensions'][1]))
        detectedWhite=np.zeros((nbins, self.myDetector.det_param['myDimensions'][0],self.myDetector.det_param['myDimensions'][1]))
        
            
        #INITIALIZING IMAGES
        incidentIntensity0=np.ones((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))*(self.exp_dict['meanShotCount']/self.exp_dict['overSampling']**2)
        incidentPhi=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        incidentDF=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        self.imageSampleBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        self.imageReferenceBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        self.imagePropagBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        self.darkFieldPropag=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
        white=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
            
        #Calculating everything for each energy of the spectrum
        for currentEnergy, flux in self.mySource.mySpectrum:
            print("Current Energy: %gkev" %currentEnergy)
            
            incidentIntensity=incidentIntensity0*flux
            if not self.exp_dict['inVacuum']:
                incidentIntensity, _, _=self.myAirVolume.setWaveRT(incidentIntensity,currentEnergy)
            
            #Take into account the detector scintillator efficiency if given in xml file
            if self.myDetector.det_param["myScintillatorMaterial"] is not None:
                for energyData, efficiency in self.myDetector.mySpectralEfficiency:
                    if energyData==currentEnergy:
                        incidentIntensity=incidentIntensity*efficiency
                
            #Passage of the incident wave through the membrane
            # print("Setting wave through membrane")
            self.IntensitySampleAfterMembrane, phiWaveSampleAfterMembrane, _=self.myMembrane.setWaveRT(incidentIntensity,currentEnergy,incidentPhi)
            
            #Propagation from membrane to object and passage through the object
            self.IntensitySampleBeforeSample,_,_=self.refraction(abs(self.IntensitySampleAfterMembrane),phiWaveSampleAfterMembrane,self.exp_dict['distMembraneToObject'],currentEnergy,self.exp_dict['magnification'])

            # print("Setting wave through sample for sample image")            
            self.IntensitySampleAfterSample,phiWaveSampleAfterSample,DarkField =self.mySampleofInterest.setWaveRT(self.IntensitySampleBeforeSample,currentEnergy,phiWaveSampleAfterMembrane,incidentDF)
            
            #Propagation to detector
            # print("Propagating waves to detector plane")
            self.imageSampleAfterRefraction,_,_=self.refraction(abs(self.IntensitySampleAfterSample),phiWaveSampleAfterSample,self.exp_dict['distObjectToDetector'],currentEnergy,self.exp_dict['magnification'], DarkField)
            self.imageReferenceAfterRefraction,_,_=self.refraction(abs(self.IntensitySampleBeforeSample),phiWaveSampleAfterMembrane,self.exp_dict['distObjectToDetector'],currentEnergy,self.exp_dict['magnification'])
            intensitySampleBeforeDetection=self.imageSampleAfterRefraction
            intensityReferenceBeforeDetection=self.imageReferenceAfterRefraction
            #Plate attenuation
            if self.myPlate is not None:
                intensitySampleBeforeDetection,_, _=self.myPlate.setWaveRT(intensitySampleBeforeDetection,currentEnergy)
                intensityReferenceBeforeDetection,_, _=self.myPlate.setWaveRT(intensityReferenceBeforeDetection,currentEnergy)
            #Combining intensities for several energies
            self.imageSampleBeforeDetection+=intensitySampleBeforeDetection
            self.imageReferenceBeforeDetection+=intensityReferenceBeforeDetection
            
            sumIntensity+=np.mean(intensityReferenceBeforeDetection)
            self.exp_dict['meanEnergy']+=currentEnergy*np.mean(intensityReferenceBeforeDetection)
            
            if pointNum==0: #We only do it for the first point
                # print("Setting wave through sample for propag and abs image")
                self.IntensityPropagAfterSample,phiWavePropagAfterSample, DarkFieldPropag=self.mySampleofInterest.setWaveRT(incidentIntensity,currentEnergy,incidentPhi, incidentDF)
                self.darkFieldPropag+=DarkFieldPropag*flux
                self.imagePropagAfterRefraction, self.Dxreal, self.Dyreal=self.refraction(abs(self.IntensityPropagAfterSample),phiWavePropagAfterSample,self.exp_dict['distObjectToDetector'],currentEnergy,self.exp_dict['magnification'], DarkFieldPropag)
                intensityPropagBeforeDetection=self.imagePropagAfterRefraction
                if self.myPlate is not None:
                    intensityPropagBeforeDetection,_,_=self.myPlate.setWaveRT(intensityPropagBeforeDetection, currentEnergy)
                    incidentIntensity,_,_=self.myPlate.setWaveRT(incidentIntensity, currentEnergy)
                white+=incidentIntensity
                self.imagePropagBeforeDetection+=intensityPropagBeforeDetection


            if currentEnergy>self.myDetector.det_param["myBinsThersholds"][ibin]-self.mySource.source_dict["myEnergySampling"]/2:
            
                effectiveSourceSize=self.mySource.source_dict["mySize"]*self.exp_dict['distObjectToDetector']/(self.exp_dict['distSourceToMembrane']+self.exp_dict['distMembraneToObject'])/self.myDetector.det_param['myPixelSize']*self.exp_dict['overSampling'] #FWHM
        
                #DETECTION IMAGES FOR ENERGY BIN
                # print("Detection sample image")
                SampleImage[ibin]=self.myDetector.detection(self.imageSampleBeforeDetection,effectiveSourceSize, self.exp_dict)
                # print("Detection reference image")
                ReferenceImage[ibin]=self.myDetector.detection(self.imageReferenceBeforeDetection,effectiveSourceSize,self.exp_dict)
                if pointNum==0:
                    self.imagePropagBeforeDetection=self.imagePropagBeforeDetection
                    # print("Detection propagation image")
                    PropagImage[ibin]=self.myDetector.detection(self.imagePropagBeforeDetection,effectiveSourceSize,self.exp_dict)
                detectedWhite[ibin]=self.myDetector.detection(white,effectiveSourceSize, self.exp_dict)
                
                self.imageSampleBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
                self.imageReferenceBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
                self.imagePropagBeforeDetection=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
                white=np.zeros((self.exp_dict['studyDimensions'][0],self.exp_dict['studyDimensions'][1]))
                
                ibin+=1
                
        self.exp_dict['meanEnergy']=self.exp_dict['meanEnergy']/sumIntensity
        print("Mean detected energy in reference image", self.exp_dict['meanEnergy'])
                
        return  SampleImage, ReferenceImage,PropagImage,detectedWhite, self.Dxreal, self.Dyreal, self.darkFieldPropag

    
        
    def saveAllParameters(self,time0,expDict):
        """
        Saves all the experimental and algorithm parameters in a txt file

        Args:
            time0 (float): time at the beginning of the calculation.
            expDict (dictionnary): dictionnary containing algorithm parameters.

        Returns:
            None.

        """
        fileName=expDict['filepath']+self.name+'_'+str(expDict['expID'])+".txt"
        print("file name: ", fileName)
        f=open(fileName,"w+")

        f.write("EXPERIMENT PARAMETERS - "+expDict['simulation_type']+" - "+str(expDict['expID']))
        
        for cle, valeur in self.exp_dict.items():
            if cle.split('_')[-1]!='unit':
                if cle+"_unit" in self.exp_dict:
                    cleUnit=cle+"_unit"
                    f.write(f'\n    {cle}: {valeur} {self.exp_dict[cleUnit]}')
                else:
                    f.write(f'\n    {cle}: {valeur}')     
            
        f.write("\n\nEntire computing time: %gs" %(time.time()-time0))   

        f.write("\n\nSource parameters:")
        f.write("\nSource name: %s" %self.mySource.myName)
        for cle, valeur in self.mySource.source_dict.items():
            if cle.split('_')[-1]!='unit':
                if cle+"_unit" in self.mySource.source_dict:
                    cleUnit=cle+"_unit"
                    f.write(f'\n    {cle}: {valeur} {self.mySource.source_dict[cleUnit]}')
                else:
                    f.write(f'\n    {cle}: {valeur}')        
        
        f.write("\n\nDetector parameters:")
        f.write("\nDetector name: %s"%self.myDetector.myName)
        for cle, valeur in self.myDetector.det_param.items():
            if cle.split('_')[-1]!='unit':
                if cle+"_unit" in self.myDetector.det_param:
                    cleUnit=cle+"_unit"
                    f.write(f'\n    {cle}: {valeur} {self.myDetector.det_param[cleUnit]}')
                else:
                    f.write(f'\n    {cle}: {valeur}')     
                    
        f.write("\n\nSample informations")
        f.write("\nSample name: %s" %self.mySampleofInterest.myName)
        f.write("\nSample type: %s" %self.mySampleType)
        f.write("\n    materials: %s" %self.mySampleofInterest.myMaterials)
        if self.mySampleofInterest.geom_parameters is not None:
            for cle, valeur in self.mySampleofInterest.geom_parameters.items():
                f.write(f'\n    {cle}: {valeur[0]} {valeur[1]}')
                
        if self.mySampleofInterest.myGeometryFunction=="openContrastPhantom":
            f.write("\nContrast Phantom geometry folder: %s" %self.mySampleofInterest.myGeometryFolder)
        
        f.write("\n\nMembrane informations:")
        f.write("\nMembrane name: %s" %self.myMembrane.myName)
        f.write("\nMembrane type: %s" %self.myMembrane.myType)
        f.write("\n    materials: %s" %self.myMembrane.myMaterials)
        f.write("\n    Membrane geometry function: %s" %self.myMembrane.myGeometryFunction)
        if self.myMembrane.geom_parameters is not None:
            for cle, valeur in self.mySampleofInterest.geom_parameters.items():
                f.write(f'\n    {cle}: {valeur[0]} {valeur[1]}')
        if self.myMembrane.myGeometryFunction=="getMembraneFromFile":
            f.write("\nMembrane geometry file: %s" %self.myMembrane.myMembraneFile)
            
        if self.myPlate is not None:
            f.write("\n\nDetectors protection Plate")
            f.write("Plate thickness: %s" %self.myPlate.myThickness)
            f.write("Plate Material: %s" %self.myPlate.myMaterials)
        
        f.close()
        
        return
    

    
        
