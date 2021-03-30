#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:14:18 2020

@author: quenot
"""

from InputOutput.pagailleIO import saveEdf,openImage
from xml.dom import minidom
import numpy as np
from Samples.createCylindre import CreateSampleCylindre
from Samples.getMembraneFromFile import getMembraneFromFile, getMembraneSegmentedFromFile
import xlrd
from matplotlib import pyplot as plt
from Samples.createSphere import CreateSampleSphere
from Samples.generateContrastPhantom import generateContrastPhantom, openContrastPhantom
import time

class Sample:
    def __init__(self):
        self.xmlSampleFileName="xmlFiles/Samples.xml"
        self.xmldocSample = minidom.parse(self.xmlSampleFileName)
        self.myName=""
        self.myType=""
        self.myMaterials=[]
        self.myGeometry=[]
        
        
    def defineCorrectValuesSample(self):
        for currentSample in self.xmldocSample.documentElement.getElementsByTagName("sample"):
            correctSample = self.getText(currentSample.getElementsByTagName("name")[0])
  
            if correctSample == self.myName:
                self.myType=(self.getText(currentSample.getElementsByTagName("myType")[0]))
                self.myMaterials=self.getText(currentSample.getElementsByTagName("myMaterials")[0])
                self.myMaterials=list(self.myMaterials.split(","))
                self.myGeometryFunction=self.getText(currentSample.getElementsByTagName("myGeometryFunction")[0])
                                
                if self.myGeometryFunction=="getMembraneFromFile" or self.myGeometryFunction=="getMembraneSegmentedFromFile":
                    self.myMembraneFile=self.getText(currentSample.getElementsByTagName("myMembraneFile")[0])
                    self.myPMMAThickness=float(self.getText(currentSample.getElementsByTagName("myPMMAThickness")[0]))
                if self.myGeometryFunction=="getMembraneSegmentedFromFile":
                    self.myMeanSphereRadius=float(self.getText(currentSample.getElementsByTagName("myMeanSphereRadius")[0]))
                    self.myNbOfLayers=int(self.getText(currentSample.getElementsByTagName("myNbOfLayers")[0]))
                    
                if self.myGeometryFunction=="get_my_thickness":
                    self.myThickness=float(self.getText(currentSample.getElementsByTagName("myThickness")[0]))
                if self.myGeometryFunction=="getVolumesFromFiles":
                    self.myGeometriesFolder=self.getText(currentSample.getElementsByTagName("myGeometriesFolder")[0])
                    self.myNumberOfSlice=int(self.getText(currentSample.getElementsByTagName("myNumberOfSlice")[0]))
                    self.myVoxelSize=float(self.getText(currentSample.getElementsByTagName("myVoxelSize")[0]))
                    self.myVolumesFiles=self.getText(currentSample.getElementsByTagName("myVolumesFiles")[0])
                    self.myVolumesFiles=list(self.myVolumesFiles.split(","))        
                    
                if self.myName=='FingerAnalytic':
                    self.myProjectionNumber=int(self.getText(currentSample.getElementsByTagName("myProjectionNumber")[0]))

                if self.myGeometryFunction=="openContrastPhantom":
                    self.myGeometryFolder=self.getText(currentSample.getElementsByTagName("myGeometryFolder")[0])

                
                return
        
        raise ValueError("Sample not found in the xml file")

    def getText(self,node):
        return node.childNodes[0].nodeValue
    
            
    def getDeltaBeta(self, sourceSpectrum):
        energyRange=[sourceSpectrum[0][0],sourceSpectrum[0][-1]]
        print("Materials :", self.myMaterials)
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
                        for energy,_ in sourceSpectrum:
                            currentCellValue=sh.cell(row,col).value
                            if energy*1000<currentCellValue:
                                delta.append((energy,0))
                                beta.append((energy,1))
                                print("No delta beta values under",currentCellValue, "eV")
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
        print("getDeltaBeta nest pas encore implemente")
            
#***********************************************************************************************************
##ANALYTICAL SAMPLE
        
class AnalyticalSample(Sample):
    def __init__(self):
        Sample.__init__(self)
        
        self.delta=[]
        self.beta=[]
        self.myGeometry=[]
        self.membranePixelSize=0
        
        
    def getMyGeometry(self,studyDimensions, studyPixelSize,oversamp, pointNum=0):
        #Returns a list of 2D arrays containing the thickness of each material of the object
        if self.myType=="sample_of_interest":
            if self.myGeometryFunction=="CreateSampleCylindre":
                self.myGeometry, self.myRadius=CreateSampleCylindre(self.myName,studyDimensions[0], studyDimensions[1], studyPixelSize)
                print("Fylon Wire Geometry")
                plt.figure()
                plt.imshow(self.myGeometry[0])
                plt.colorbar()
                plt.show()
                return                
            if self.myGeometryFunction=="generateContrastPhantom":
                self.myGeometry=generateContrastPhantom(studyDimensions[0],studyDimensions[1],studyPixelSize, angle=90)
                return
            if self.myGeometryFunction=="openContrastPhantom":
                self.myGeometry=openContrastPhantom(self.myGeometryFolder,studyDimensions[0],studyDimensions[1],studyPixelSize,oversamp, angle=90)
                return
        
        if self.myType=="membrane":
            if self.myGeometryFunction=="getMembraneFromFile":
                self.myGeometry.append(getMembraneFromFile(self.myMembraneFile,studyDimensions,studyPixelSize,oversamp,pointNum))
                self.myGeometry.append(np.ones((studyDimensions[0], studyDimensions[1]))*self.myPMMAThickness*1e-6)
                return
            if self.myGeometryFunction=="getMembraneSegmentedFromFile":
                self.myGeometry.append(getMembraneSegmentedFromFile(self,studyDimensions[0],studyDimensions[1],studyPixelSize,oversamp,pointNum))
                self.myGeometry.append(np.ones((studyDimensions[0], studyDimensions[1]))*self.myPMMAThickness*1e-6)
                return
        
        if self.myGeometryFunction=="get_my_thickness":
            self.myGeometry=[]
            self.myGeometry.append(np.ones((studyDimensions[0], studyDimensions[1]))*self.myThickness*1e-6)
            return
        
        raise ValueError("Could not define sample geometry")
        
        
    def setWave(self,incidentWave, energy):
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
                    print(energyData, betaEn)
            disturbedWave=np.exp((-1j*k*delta[imat]-k*beta[imat])*self.myGeometry[imat])*disturbedWave
            
        print('Finished set wave through', self.myName)
        return disturbedWave
    
            
    def setWaveRT(self,incidentIntensity,incidentphi, energy, distFromSource):
        #Returns the intensity map after absorption in the object and the phase shift (specific to ray tracing)
        
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
#
            geometry=self.myGeometry[imat]
            disturbedIntensity=np.exp(-2*k*beta[imat]*self.myGeometry[imat])*disturbedIntensity
            disturbedPhi=disturbedPhi-k*delta[imat]*self.myGeometry[imat]
            
        print('Finished setting wave through', self.myName)
        return disturbedIntensity, disturbedPhi
        
        