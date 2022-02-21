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
from numpy.fft import fftshift, ifftshift, fft2, ifft2
from numpy import pi
from refractionFileNumba import fastRefraction
from getk import getk
import multiprocess as mp
from minOversampling import minOversampling

class Experiment:
    def __init__(self, exp_dict):
        """Experiment class constructor.
            
        Args:
            exp_dict (dict): dictionary of simulation algorithm parameters
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
        self.simulation_type=exp_dict['simulation_type']
        self.meanShotCount=0
        self.meanEnergy=0
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
        
        if self.myDetector.myScintillatorMaterial is not None:
            self.myDetector.getBeta(self.mySource.mySpectrum)

        min_oversampling = max([minOversampling(self.distObjectToDetector, self.distSourceToMembrane+self.distMembraneToObject, energy, self.myDetector.myPixelSize) for energy,_ in self.mySource.mySpectrum]) \
            if self.simulation_type == 'Fresnel' else 2
        # if self.sampling is None: #TODO get this bit working so that None gives the minimum oversampling
        #     self.sampling = min_oversampling
        if self.sampling < min_oversampling:
            raise ValueError(f'Current sampling ({self.sampling}) is insufficient. The minimum sampling required is: {min_oversampling}')
        
        print(self.mySource)
        print(self.myDetector)
        print(self.mySampleofInterest)
        print(self.myMembrane)
        print(f'Study:\n'
        f' Magnification: {self.magnification}\n'
        f' Study Dimensions: {self.studyDimensions} pixels\n'
        f' Study Pixel Size: {self.studyPixelSize} um\n'
        f' Oversampling: {self.sampling}\n'
        f' Minimum Oversampling: {min_oversampling}\n')

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
        """
        Calculates the study dimensions considereing the geometry of the set up, the field of view and the sample pixels oversampling

        Returns:
            None.

        """
        self.precision=(self.myDetector.myPixelSize/self.sampling/self.distObjectToDetector)
        self.studyDimensions=self.myDetector.myDimensions*int(self.sampling)
        self.studyPixelSize=self.myDetector.myPixelSize/self.sampling/self.magnification
        
    
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
        if propagationDistance == 0:
            return waveToPropagate
        
        #Propagateur de Fresnel
        k = getk(Energy*1000)
        
        Nx, Ny = self.studyDimensions     
        u, v = np.meshgrid(np.arange(0, Nx), np.arange(0, Ny))
        u = (u - (Nx / 2))
        v = (v - (Ny / 2))
        u_m = u *2*pi / (self.studyDimensions[0]*self.studyPixelSize*1e-6)
        v_m = v *2*pi / (self.studyDimensions[1]*self.studyPixelSize*1e-6)
        uv_sqr=  np.transpose(u_m ** 2 + v_m ** 2)  # ie (u2+v2)
        
        waveAfterPropagation=np.exp(1j*k*propagationDistance/magnification)*ifft2(ifftshift(np.exp(-1j*propagationDistance*(uv_sqr)/(2*k*magnification))*fftshift(fft2(waveToPropagate))))
        
        return waveAfterPropagation
    
    
    def refraction(self, intensityRefracted, phi, propagationDistance, Energy, magnification):
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
        intensityRefracted2, Dx, Dy = fastRefraction(intensityRefracted, phi, propagationDistance, Energy, magnification, self.studyPixelSize)
        return intensityRefracted2,Dx, Dy    
    
    def computeSampleAndReferenceImages(self, pointNum):
        """
        Compute intensity changes on the path of the previously difined experiment 
        to create all the images of the SBI experiment with Fresnel propagator


        Returns:
            SampleImage (2d numpy array): sample image simulated with sample and membrane.
            ReferenceImage (2d numpy array): reference image simulated with membrane.
            PropagImage (2d numpy array): propagation image simulated with only sample.
            detectedWhite (2d numpy array): white image without membrane and sample.

        """
        
        #INITIALIZING PARAMETERS
        if pointNum==0:
            if any(elem<self.mySource.mySpectrum[0][0] for elem in self.myDetector.myBinsThresholds) or any(elem>self.mySource.mySpectrum[-1][0] for elem in self.myDetector.myBinsThresholds):
                raise Exception(f'At least one of your detector bin threshold is outside your source spectrum. \nYour source spectrum ranges from {self.mySource.mySpectrum[0][0]} to {self.mySource.mySpectrum[-1][0]}')
            # self.myDetector.myBinsThresholds.insert(0,self.mySource.mySpectrum[0][0])
            # if self.myDetector.myBinsThresholds[-1] != self.mySource.mySpectrum[-1][0]:
            #     self.myDetector.myBinsThresholds.append(self.mySource.mySpectrum[-1][0])
    
        effectiveSourceSize=self.mySource.mySize*self.distObjectToDetector/(self.distSourceToMembrane+self.distMembraneToObject)/self.myDetector.myPixelSize*self.sampling #FWHM
        
        #Defining total flux for normalizing spectrum
        energies, fluxes = list(zip(*[(currentEnergy, flux) for currentEnergy, flux in self.mySource.mySpectrum]))

        #####################################################################################
        global _parallel_propagate
        #Calculating everything for each energy of the spectrum
        def _parallel_propagate(currentEnergy, flux):
            print("\nCurrent Energy:", currentEnergy)
            #Taking into account source window and air attenuation of intensity
            incidentIntensity=np.ones([*self.studyDimensions])*self.meanShotCount*flux/self.mySource.totalFlux()
            incidentIntensity, _=self.myAirVolume.setWaveRT(incidentIntensity,1, currentEnergy)
            
            #Take into account the detector scintillator efficiency if given in xml file
            if self.myDetector.myScintillatorMaterial is not None:
                beta = [betaEn for energyData, betaEn  in self.myDetector.beta if energyData==currentEnergy][0]
                k=getk(currentEnergy*1e3)
                detectedSpectrum=1-np.exp(-2*k*self.myDetector.myScintillatorThickness*1e-6*beta)
                print("Scintillator efficiency for current energy:", detectedSpectrum)
                incidentIntensity=incidentIntensity*detectedSpectrum
            incidentWave=np.sqrt(incidentIntensity)
            
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
            
            
            Intensity = np.mean(intensityReferenceBeforeDetection)
            self.meanEnergy+=currentEnergy*Intensity
            
            if pointNum==0: #We only do it for the first point
                print("Setting wave through sample for propag and abs image")
                self.wavePropagAfterSample=self.mySampleofInterest.setWave(incidentWave,currentEnergy)
                self.wavePropagBeforeDetection=self.wavePropagation(self.wavePropagAfterSample,self.distObjectToDetector,currentEnergy,self.magnification)
                intensityPropagBeforeDetection=abs(self.wavePropagBeforeDetection)**2
                if self.myPlaque is not None:
                    intensityPropagBeforeDetection,_=self.myPlaque.setWaveRT(intensityPropagBeforeDetection,1, currentEnergy)
                    #######
                incidentIntensityWhite=incidentWave**2
                if self.myPlaque is not None:
                    incidentIntensityWhite,_=self.myPlaque.setWaveRT(incidentWave**2,1, currentEnergy)
            # XXX Not sure if this is necessary???
            else:
                intensityPropagBeforeDetection = incidentIntensityWhite = np.zeros([*self.studyDimensions])

            return intensitySampleBeforeDetection, intensityReferenceBeforeDetection, intensityPropagBeforeDetection, incidentIntensityWhite

        splits = np.searchsorted(energies, self.myDetector.myBinsThresholds, side='right')
        binned_energies = np.split(energies, splits)
        binned_fluxes = np.split(fluxes, splits)
        binned_energies = [i for i in binned_energies if i.size]
        binned_fluxes = [i for i in binned_fluxes if i.size]
        self.imageSampleBeforeDetection = []
        self.imageReferenceBeforeDetection = []
        self.imagePropagBeforeDetection = []
        self.white = []
        # TODO only do it with pool if multiprocessing, so only have one experiment file
        # TODO pass the cpu count as an argument
        pool = mp.Pool(processes = mp.cpu_count())
        for bin_en, bin_flu in zip(binned_energies, binned_fluxes):
            binlength = len(bin_en)
            isbd = np.zeros((binlength, *self.studyDimensions))
            irbd = np.zeros((binlength, *self.studyDimensions))
            ipbd = np.zeros((binlength, *self.studyDimensions))
            iw = np.zeros((binlength, *self.studyDimensions))

            sum_isbd = np.zeros(self.studyDimensions)
            sum_irbd = np.zeros(self.studyDimensions)
            sum_ipbd = np.zeros(self.studyDimensions)
            sum_iw = np.zeros(self.studyDimensions)

            isbd, irbd, ipbd, iw = np.array(list(zip(*np.array(pool.starmap(_parallel_propagate, zip(bin_en, bin_flu))))))
            sum_isbd = isbd.sum(0)
            sum_irbd = irbd.sum(0)
            sum_ipbd = ipbd.sum(0)
            sum_iw = iw.sum(0)
            self.imageSampleBeforeDetection.append(sum_isbd)
            self.imageReferenceBeforeDetection.append(sum_irbd)
            self.imagePropagBeforeDetection.append(sum_ipbd)
            self.white.append(sum_iw)
        
        SampleImage = []
        ReferenceImage = []
        PropagImage = []
        detectedWhite = []

        for sample, reference, propagation, white in zip(self.imageSampleBeforeDetection,
            self.imageReferenceBeforeDetection, self.imagePropagBeforeDetection, self.white):

            print("Detection sample image")
            SampleImage.append(self.myDetector.detection(sample, effectiveSourceSize))
            print("Detection reference image")
            ReferenceImage.append(self.myDetector.detection(reference,effectiveSourceSize))
            if pointNum==0:
                print("Detection propagation image")
                PropagImage.append(self.myDetector.detection(propagation,effectiveSourceSize))
            detectedWhite.append(self.myDetector.detection(white, effectiveSourceSize))

        return SampleImage, ReferenceImage, PropagImage, detectedWhite

        
    def computeSampleAndReferenceImagesRT(self, pointNum):
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
        if pointNum==0:
            if any(elem<self.mySource.mySpectrum[0][0] for elem in self.myDetector.myBinsThresholds) or any(elem>self.mySource.mySpectrum[-1][0] for elem in self.myDetector.myBinsThresholds):
                raise Exception(f'At least one of your detector bin Threshold is outside your source spectrum. \nYour source spectrum ranges from {self.mySource.mySpectrum[0][0]} to {self.mySource.mySpectrum[-1][0]}')
            # self.myDetector.myBinsThresholds.insert(0,self.mySource.mySpectrum[0][0])
            self.myDetector.myBinsThresholds.append(self.mySource.mySpectrum[-1][0])

        effectiveSourceSize=self.mySource.mySize*self.distObjectToDetector/(self.distSourceToMembrane+self.distMembraneToObject)/self.myDetector.myPixelSize*self.sampling #FWHM

        energies, fluxes = list(zip(*[(currentEnergy, flux) for currentEnergy, flux in self.mySource.mySpectrum]))
        incidentPhi = np.zeros([*self.studyDimensions])

        global _parallel_propagateRT

        def _parallel_propagateRT(currentEnergy, flux):

            #Calculating everything for each energy of the spectrum
            print("\nCurrent Energy: %gkev" %currentEnergy)
            
            incidentIntensity=np.ones([*self.studyDimensions])*self.meanShotCount*flux/self.mySource.totalFlux()
            incidentIntensityWhite,_ =self.myAirVolume.setWaveRT(incidentIntensity,1, currentEnergy)
            
            #Take into account the detector scintillator efficiency if given in xml file
            if self.myDetector.myScintillatorMaterial is not None:
                beta = [betaEn for energyData, betaEn  in self.myDetector.beta if energyData==currentEnergy][0]
                k=getk(currentEnergy*1e3)
                detectedSpectrum=1-np.exp(-2*k*self.myDetector.myScintillatorThickness*1e-6*beta)
                print("Scintillator efficiency for current energy:", detectedSpectrum)
                incidentIntensity=incidentIntensity*detectedSpectrum
                
            #Passage of the incident wave through the membrane
            print("Setting wave through membrane")
            self.IntensitySampleAfterMembrane, phiWaveSampleAfterMembrane=self.myMembrane.setWaveRT(incidentIntensity,incidentPhi,currentEnergy)
            
            #Propagation from membrane to object and passage through the object
            self.IntensitySampleBeforeSample,_,_=self.refraction(abs(self.IntensitySampleAfterMembrane),phiWaveSampleAfterMembrane,self.distMembraneToObject,currentEnergy,self.magnification)

            print("Setting wave through sample for sample image")            
            self.IntensitySampleAfterSample,phiWaveSampleAfterSample=self.mySampleofInterest.setWaveRT(self.IntensitySampleBeforeSample,phiWaveSampleAfterMembrane,currentEnergy)
            
            #Propagation to detector
            print("Propagating waves to detector plane")
            intensitySampleBeforeDetection,_,_=self.refraction(abs(self.IntensitySampleAfterSample),phiWaveSampleAfterSample,self.distObjectToDetector,currentEnergy,self.magnification)
            intensityReferenceBeforeDetection,_,_=self.refraction(abs(self.IntensitySampleBeforeSample),phiWaveSampleAfterMembrane,self.distObjectToDetector,currentEnergy,self.magnification)
            #Plaque attenuation
            if self.myPlaque is not None:
                intensitySampleBeforeDetection,_=self.myPlaque.setWaveRT(intensitySampleBeforeDetection,1, currentEnergy)
                intensityReferenceBeforeDetection,_=self.myPlaque.setWaveRT(intensityReferenceBeforeDetection,1, currentEnergy)
            
            if pointNum==0: #We only do it for the first point
                print("Setting wave through sample for propag and abs image")
                self.IntensityPropagAfterSample,phiWavePropagAfterSample=self.mySampleofInterest.setWaveRT(incidentIntensity,incidentPhi,currentEnergy)
                self.imagePropagAfterRefraction, self.Dxreal, self.Dyreal=self.refraction(abs(self.IntensityPropagAfterSample),phiWavePropagAfterSample,self.distObjectToDetector,currentEnergy,self.magnification)
                intensityPropagBeforeDetection=self.imagePropagAfterRefraction
                if self.myPlaque is not None:
                    intensityPropagBeforeDetection,_=self.myPlaque.setWaveRT(intensityPropagBeforeDetection,1, currentEnergy)
                    incidentIntensityWhite,_=self.myPlaque.setWaveRT(incidentIntensity,1, currentEnergy)
            #  XXX Again, not sure if this is necessary???
            else:
                intensityPropagBeforeDetection = incidentIntensityWhite = np.zeros([*self.studyDimensions])

            return intensitySampleBeforeDetection, intensityReferenceBeforeDetection, intensityPropagBeforeDetection, incidentIntensityWhite

        splits = np.searchsorted(energies, self.myDetector.myBinsThresholds, side='right')
        binned_energies = np.split(energies, splits)
        binned_fluxes = np.split(fluxes, splits)
        binned_energies = [i for i in binned_energies if i.size]
        binned_fluxes = [i for i in binned_fluxes if i.size]
        self.imageSampleBeforeDetection = []
        self.imageReferenceBeforeDetection = []
        self.imagePropagBeforeDetection = []
        self.white = []
        # TODO only do it with pool if multiprocessing, so only have one experiment file
        # TODO pass the cpu count as an argument
        pool = mp.Pool(processes = mp.cpu_count())
        for bin_en, bin_flu in zip(binned_energies, binned_fluxes):
            binlength = len(bin_en)
            isbd = np.zeros((binlength, *self.studyDimensions))
            irbd = np.zeros((binlength, *self.studyDimensions))
            ipbd = np.zeros((binlength, *self.studyDimensions))
            iw = np.zeros((binlength, *self.studyDimensions))

            sum_isbd = np.zeros(self.studyDimensions)
            sum_irbd = np.zeros(self.studyDimensions)
            sum_ipbd = np.zeros(self.studyDimensions)
            sum_iw = np.zeros(self.studyDimensions)

            isbd, irbd, ipbd, iw = np.array(list(zip(*np.array(pool.starmap(_parallel_propagateRT, zip(bin_en, bin_flu))))))
            sum_isbd = isbd.sum(0)
            sum_irbd = irbd.sum(0)
            sum_ipbd = ipbd.sum(0)
            sum_iw = iw.sum(0)
            self.imageSampleBeforeDetection.append(sum_isbd)
            self.imageReferenceBeforeDetection.append(sum_irbd)
            self.imagePropagBeforeDetection.append(sum_ipbd)
            self.white.append(sum_iw)
        
        SampleImage = []
        ReferenceImage = []
        PropagImage = []
        detectedWhite = []

        for sample, reference, propagation, white in zip(self.imageSampleBeforeDetection,
            self.imageReferenceBeforeDetection, self.imagePropagBeforeDetection, self.white):

            print("Detection sample image")
            SampleImage.append(self.myDetector.detection(sample, effectiveSourceSize))
            print("Detection reference image")
            ReferenceImage.append(self.myDetector.detection(reference,effectiveSourceSize))
            if pointNum==0:
                print("Detection propagation image")
                PropagImage.append(self.myDetector.detection(propagation,effectiveSourceSize))
            detectedWhite.append(self.myDetector.detection(white, effectiveSourceSize))

        return SampleImage, ReferenceImage, PropagImage, detectedWhite
    
        
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

        f.write("\n\nDistances:")
        f.write("\nDistance source to membrane: %gm" %self.distSourceToMembrane)
        f.write("\nDistance membrane to sample: %gm" %self.distMembraneToObject)
        f.write("\nDistance sample to detector: %gm" %self.distObjectToDetector)

        f.write("\n\nSource parameters:")
        f.write("\nSource name: %s" %self.mySource.myName)
        f.write("\nSource type: %s" %self.mySource.myType)
        f.write("\nSource size: %gum" %self.mySource.mySize)
        if self.mySource.myType=="Monochromatic":
            f.write("\nSource energy: %gkev" %(self.mySource.mySpectrum[0][0]))
        if self.mySource.myType=="Polychromatic":
            f.write("\nSource voltage: %gkVp" %self.mySource.myVoltage)
            f.write("\nSource spectrum energy sampling: %gkeV" %self.mySource.myEnergySampling)
            f.write("\nMean energy detected in the reference image: %gkeV" %self.meanEnergy)
        
        f.write("\n\nDetector parameters:")
        f.write("\nDetector name: %s"%self.myDetector.myName)
        f.write("\nDetector dimensions:"+str(self.myDetector.myDimensions[0]-expDict['margin'])+"x"+str(self.myDetector.myDimensions[1]-expDict['margin'])+"pix")
        f.write("\nDetector pixels size: %gum" %self.myDetector.myPixelSize)
        f.write("\nDetector PSF: %gpix" %self.myDetector.myPSF)
        if self.myDetector.myEnergyLimit is not None:
            f.write("\nDetector max energy detected: %gkeV" %self.myDetector.myEnergyLimit)
        if self.myDetector.myBinsThresholds is not None:
            f.write(f'\nDetector bins Thresholds: {self.myDetector.myBinsThresholds}keV')
        if self.myDetector.myScintillatorMaterial is not None:
            f.write(f"  Scintillator {self.myDetector.myScintillatorMaterial} of {self.myDetector.myScintillatorThickness}um")
        
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
        f.write("\nMean shot count: %g" %(self.meanShotCount*self.sampling**2))        
        f.write("\nNumber of points: %g" %(expDict['nbExpPoints']))    
        f.write("\nEntire computing time: %gs" %(time.time()-time0))   
        
        f.close()
        
        return
    