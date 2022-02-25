#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:14:18 2020

@author: quenot
"""

from InputOutput.pagailleIO import saveEdf,openImage
from xml.dom import minidom
import numpy as np
from Samples.getMembraneFromFile import getMembraneFromFile, getMembraneSegmentedFromFile
import xlrd
from matplotlib import pyplot as plt
from Samples.generateContrastPhantom import generateContrastPhantom
from Samples.createSampGeom import CreateSampleSpheresInParallelepiped, CreateSampleSpheresInCylinder, CreateSampleCylindre, CreateSampleSphere, CreateYourSampleGeometry, loadSampleGeometryFromImages
import time
from numba import jit 

class Sample:
    def __init__(self):
        self.xmlSampleFileName="xmlFiles/Samples.xml"
        self.xmldocSample = minidom.parse(self.xmlSampleFileName)
        self.myName=""
        self.myType=""
        self.myMaterials=[]
        self.myGeometry=[]
        self.geom_parameters=None
        
        
    def defineCorrectValuesSample(self):
        """
        gets all sample parameters from xml file

        Raises:
            ValueError: Sample not found in the xml file.

        Returns:
            None.

        """
        for currentSample in self.xmldocSample.documentElement.getElementsByTagName("sample"):
            correctSample = self.getText(currentSample.getElementsByTagName("name")[0])
  
            if correctSample == self.myName:
                self.myType=(self.getText(currentSample.getElementsByTagName("myType")[0]))
                self.myMaterials=self.getText(currentSample.getElementsByTagName("myMaterials")[0])
                self.myMaterials=list(self.myMaterials.split(","))
                self.myGeometryFunction=self.getText(currentSample.getElementsByTagName("myGeometryFunction")[0])
                                
                if self.myGeometryFunction=="getMembraneFromFile":
                    self.myPMMAThickness=float(self.getText(currentSample.getElementsByTagName("myPMMAThickness")[0]))
                if self.myGeometryFunction=="getMembraneSegmentedFromFile":
                    self.myMeanSphereRadius=float(self.getText(currentSample.getElementsByTagName("myMeanSphereRadius")[0]))
                    self.myNbOfLayers=int(self.getText(currentSample.getElementsByTagName("myNbOfLayers")[0]))
                    self.myPMMAThickness=float(self.getText(currentSample.getElementsByTagName("myPMMAThickness")[0]))
                
                if self.myGeometryFunction=="get_my_thickness":
                    if self.myName!="air_volume":
                        self.myThickness=float(self.getText(currentSample.getElementsByTagName("myThickness")[0]))
                if self.myGeometryFunction=="getSampleFromFile":
                    self.mySampleFile=self.getText(currentSample.getElementsByTagName("mySampleFile")[0])
                if self.myGeometryFunction=="getVolumesFromFiles":
                    self.myGeometriesFolder=self.getText(currentSample.getElementsByTagName("myGeometriesFolder")[0])
                    self.myNumberOfSlice=int(self.getText(currentSample.getElementsByTagName("myNumberOfSlice")[0]))
                    self.myVoxelSize=float(self.getText(currentSample.getElementsByTagName("myVoxelSize")[0]))
                    self.myVolumesFiles=self.getText(currentSample.getElementsByTagName("myVolumesFiles")[0])
                    self.myVolumesFiles=list(self.myVolumesFiles.split(","))        
                
                if self.myGeometryFunction=="loadSampleGeometryFromImages":
                    self.myGeometryFolder=self.getText(currentSample.getElementsByTagName("myGeometryFolder")[0])

                return
        print(self.myName)
        raise ValueError("Sample not found in the xml file")

    def getText(self,node):
        return node.childNodes[0].nodeValue
    
            
    def getDeltaBeta(self, sourceSpectrum):
        """
        gets each materials delta and beta parameters for every energy of the spectrum

        Args:
            sourceSpectrum (list of tuples): spectrum of the source containing energy (in keV) and weights.

        Returns:
            None.

        """
        energyRange=[sourceSpectrum[0][0],sourceSpectrum[0][-1]]
        pathTablesDeltaBeta ='Samples/DeltaBeta/TablesDeltaBeta.xls'
        deltaBetaDoc=xlrd.open_workbook(pathTablesDeltaBeta)
        
        i=0
        for sh in xlrd.open_workbook(pathTablesDeltaBeta).sheets():
            for imat in range(len(self.myMaterials)):
                delta=[]
                beta=[]
                for col in range(sh.ncols):
                    row=0
                    myCell = sh.cell(row, col)
#                    print("\n\n Sample materials : ",self.myMaterials[imat])
                    if myCell.value == self.myMaterials[imat]:
                        row=row+3
                        lowEnergyAbsent=False
                        for energy,_ in sourceSpectrum:
                            currentCellValue=sh.cell(row,col).value
                            if energy*1000<currentCellValue: #E is in keV in the spectrum and in eV in the tables delta-beta file
                                delta.append((energy,0))
                                beta.append((energy,1))
                                if not lowEnergyAbsent:
                                    print(f'/!\ No delta beta values under {currentCellValue/1000} keV for {self.myMaterials[imat]}')
                                    lowEnergyAbsent=True
                                continue
                            nextCellValue=sh.cell(row+1,col).value
                            while nextCellValue<energy*1e3: #find in which interval the energy is in the delta beta tables
                                row+=1
                                currentCellValue=sh.cell(row,col).value
                                nextCellValue=sh.cell(row+1,col).value
                            #Linear interpollation between those values
                            step=nextCellValue-currentCellValue
                            currentCellDelta=float(sh.cell(row,col+1).value)
                            nextCellDelta=float(sh.cell(row+1,col+1).value)
                            currentCellBeta=float(sh.cell(row,col+2).value)
                            nextCellBeta=float(sh.cell(row+1,col+2).value)
                            deltaInterp=abs(nextCellValue-energy*1e3)/step*currentCellDelta+abs(currentCellValue-energy*1e3)/step*nextCellDelta
                            betaInterp=abs(nextCellValue-energy*1e3)/step*currentCellBeta+abs(currentCellValue-energy*1e3)/step*nextCellBeta
                            delta.append((energy,deltaInterp))
                            beta.append((energy,betaInterp))
                               
                            
                        self.beta.append(beta)
                        self.delta.append(delta)
                        
            a, b, c=np.shape(self.delta)
            if a<len(self.myMaterials):
                raise ValueError("One or more materials have not been found in delta beta tables")
            return
            
#***********************************************************************************************************
##ANALYTICAL SAMPLE
        
class AnalyticalSample(Sample):
    def __init__(self):
        Sample.__init__(self)
        self.delta=[]
        self.beta=[]        
        
    def getMyGeometry(self,studyDimensions, studyPixelSize,oversamp, pointNum=0, number_of_positions=0):
        """
        Sets sample geometry maps for each of their materials geometry[material, x, y]

        Args:
            studyDimensions (tuple): (Dimx, Dimy).
            studyPixelSize (float): voxel size in the considered plane in um (taking into account the overoverSampling).
            oversamp (int): overoverSampling compared to detector pixels.
            pointNum (int, optional): number of the calculated point (only useful to move the membrane between each point). Defaults to 0.

        Raises:
            ValueError: Could not define sample geometry means that the case was not found.

        Returns:
            None.

        """
        #Returns a list of 2D arrays containing the thickness of each material of the object
        if self.myType=="sample_of_interest":
            if self.myGeometryFunction=="getSampleFromFile":
                self.myGeometry = np.load(self.mySampleFile)
                plt.figure()
                plt.imshow(self.myGeometry[0])
                plt.title()
                plt.colorbar("Sample Geometry")
                plt.show()
                return   
            if self.myGeometryFunction=="CreateSampleCylindre":
                self.myGeometry,self.geom_parameters=CreateSampleCylindre(self.myName,studyDimensions[0], studyDimensions[1], studyPixelSize)
                plt.figure()
                plt.imshow(self.myGeometry[0])
                plt.title("Cylinder Geometry (m)")
                plt.colorbar()
                plt.show()
                return                
            if self.myGeometryFunction=="CreateYourSampleGeometry":
                self.myGeometry,self.geom_parameters=CreateYourSampleGeometry(self.myName,studyDimensions[0], studyDimensions[1], studyPixelSize)
                plt.figure()
                plt.imshow(self.myGeometry[0])
                plt.title("Geometry (m)")
                plt.colorbar()
                plt.show()
                return                
            if self.myGeometryFunction=="CreateSampleSpheresInCylinder":
                self.myGeometry,self.geom_parameters =CreateSampleSpheresInCylinder(self.myName,studyDimensions[0], studyDimensions[1], studyPixelSize)
                return              
            if self.myGeometryFunction=="CreateSampleSpheresInParallelepiped":
                self.myGeometry,self.geom_parameters =CreateSampleSpheresInParallelepiped(self.myName,studyDimensions[0], studyDimensions[1], studyPixelSize)
                return              
            if self.myGeometryFunction=="CreateSampleSphere":
                self.myGeometry,self.geom_parameters =CreateSampleSphere(self.myName,studyDimensions[0], studyDimensions[1], studyPixelSize)
                plt.figure()
                plt.imshow(self.myGeometry[0])
                plt.title("Sphere geometry (m)")
                plt.colorbar()
                plt.show()
                return       
            if self.myGeometryFunction=="generateContrastPhantom":
                self.myGeometry,self.geom_parameters=generateContrastPhantom(studyDimensions[0],studyDimensions[1],studyPixelSize, angle=30)
                self.myGeometry=np.array(self.myGeometry)
                return
            if self.myGeometryFunction=="loadSampleGeometryFromImages":
                self.myGeometry,self.geom_parameters=loadSampleGeometryFromImages(self.myGeometryFolder,studyDimensions[0],studyDimensions[1],studyPixelSize)
                self.myGeometry=np.array(self.myGeometry)
                return
        
        if self.myType=="membrane":
            if self.myGeometryFunction=="getMembraneFromFile":
                self.myGeometry,self.geom_parameters=getMembraneFromFile(self.myMembraneFile,studyDimensions,pointNum, self.myPMMAThickness)
                self.myGeometry=np.array(self.myGeometry)
                return
            if self.myGeometryFunction=="getMembraneSegmentedFromFile":
                self.myGeometry,self.geom_parameters=getMembraneSegmentedFromFile(self,studyDimensions[0],studyDimensions[1],studyPixelSize,pointNum, self.myPMMAThickness)
                self.myGeometry=np.array(self.myGeometry)
                return
        
        if self.myGeometryFunction=="get_my_thickness":
            self.myGeometry=[]
            self.myGeometry.append(np.ones((studyDimensions[0], studyDimensions[1]))*self.myThickness*1e-6)
            self.myGeometry=np.array(self.myGeometry)
            return
        
        raise ValueError("Could not define sample geometry")
        
        
    def setWave(self,incidentWave, energy):
        """
        sets waves through the sample

        Args:
            incidentWave (2d numpy array): wave before the sample.
            energy (float): considered energy in keV.

        Raises:
            Exception: Sample Geometry has the wrong nb of dim [material, x, y].

        Returns:
            disturbedWave (2d numpy array): wave after the sample.

        """
        if self.myGeometry.ndim != 3:
            raise Exception("Sample Geometry has the wrong nb of dim [material, x, y]")
        k=2*np.pi*energy*1000*1.6e-19/(6.626e-34*2.998e8)
        delta=np.zeros(len(self.myMaterials))
        beta=np.zeros(len(self.myMaterials))
                    
        disturbedWave=incidentWave
        for imat in range(len(self.myMaterials)):
            deltaList=self.delta[imat]
            for energyData, deltaEn in deltaList:
                if energyData==energy:
                    delta[imat]=deltaEn
            for energyData, betaEn in self.beta[imat]:
                if energyData==energy:
                    beta[imat]=betaEn
                    # print(energyData, betaEn)
            disturbedWave=np.exp((-1j*k*delta[imat]-k*beta[imat])*self.myGeometry[imat])*disturbedWave
            
        # print('Finished set wave through', self.myName)
        return disturbedWave
    
            
    def setWaveRT(self,incidentIntensity,energy,incidentphi=0,incidentDf=0):
        """
        sets intensity and phase through the sample

        Args:
            incidentIntensity (2d numpy array): DESCRIPTION.
            incidentphi (2d numpy array): DESCRIPTION.
            energy (float): current energy in keV.

        Returns:
            disturbedIntensity (2d numpy array): intensity after attenuation of the sample.
            disturbedPhi (2d numpy array): phase shift due to the sample.

        """        
        newDf=0
        k=2*np.pi*energy*1000*1.6e-19/(6.626e-34*2.998e8)

        #Get delta and beta of the materials for the considered energy
        delta=np.zeros(len(self.myMaterials))
        beta=np.zeros(len(self.myMaterials))
        
        geometryBefore=self.myGeometry
                
        disturbedIntensity=incidentIntensity
        disturbedPhi=incidentphi
        for imat in range(len(self.myMaterials)):
            deltaList=self.delta[imat]
            
            #retrieve delta and beta from list of the material considered
            for energyData, deltaEn in deltaList:
                if energyData==energy:
                    delta[imat]=deltaEn
            for energyData, betaEn in self.beta[imat]:
                if energyData==energy:
                    beta[imat]=betaEn
            geometry=self.myGeometry[imat]
#
            if self.myType=="sample_of_interest":
                
                if self.myMaterials[imat]=="Lung":
                    alveoliRadius=47 #um
                    sphereVolumeFraction=0.5
                    deltaDelta=delta[imat]
                    geom=self.myGeometry[imat]*1e6 #um volume actually filled with wolume against air
                    NsphereVol=sphereVolumeFraction*3/4/np.pi/(alveoliRadius**3)
                    Nsphere=NsphereVol**(1/3)*geom#3/4*geom/alveoliRadius
                    newDf=2*deltaDelta*(Nsphere)**(1/2)*np.sqrt(np.log(2/deltaDelta)+1)
                    geometry=geometry*sphereVolumeFraction
                    # newDf[newDf<0.01]=0
                
                if self.myName=='cylinder_beeds':
                    alveoliRadius=15 #um
                    sphereVolumeFraction=0.6
                    deltaDelta=delta[imat]
                    geom=self.myGeometry[imat]*1e6 #um volume actually filled with wolume against air
                    NsphereVol=sphereVolumeFraction*3/4/np.pi/(alveoliRadius**3)
                    Nsphere=NsphereVol**(1/3)*geom#3/4*geom/alveoliRadius
                    newDf=2*deltaDelta*(Nsphere)**(1/2)*np.sqrt(np.log(2/deltaDelta)+1)
                    geometry=geometry*sphereVolumeFraction
                    # newDf[newDf<0.01]=0
                
            
            disturbedIntensity=np.exp(-2*k*beta[imat]*geometry)*disturbedIntensity
            disturbedPhi=disturbedPhi-k*delta[imat]*geometry
                                    
        # print('Finished setting wave through', self.myName)
        return disturbedIntensity, disturbedPhi, newDf
        
        
