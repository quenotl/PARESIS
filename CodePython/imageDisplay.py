#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 11:17:27 2020

@author: quenot
"""
from matplotlib import pyplot as plt

def imageDisplay(experiment,SampleImage, ReferenceImage,PropagImage, dxOF, dyOF, phiLarkin, SubImage ):
    
        
    ## Useful displays
    SampleImafterSample=abs(experiment.waveSampleAfterSample**2)#.imag#(abs(waveSampleAfterSample))**2
    sampleImBeforeDetection=experiment.imageSampleBeforeDetection
    ReferenceImafterMembrane=abs(experiment.waveSampleAfterMembrane**2)#.imag#(abs(waveSampleAfterSample))**2
    referenceImBeforeDetection=experiment.imageReferenceBeforeDetection

    
    print('\Membrane Thickness (m)')
    plt.figure()
    plt.imshow(experiment.myMembrane.myGeometry[0])
    plt.colorbar()
    plt.show()
    
    print('\ReferenceImafterMembrane')
    plt.figure()
    plt.imshow(ReferenceImafterMembrane)
    plt.colorbar()
    plt.show()
        
    print('\Im sample after object')
    plt.figure()
    plt.imshow(SampleImafterSample)
    plt.colorbar()
    plt.show()

#    print('\nIm sample before detection')
#    plt.figure()
#    plt.imshow(sampleImBeforeDetection[0])
#    plt.colorbar()
#    plt.show()
    
    print('\nDetected Sample Image')
    plt.figure()
    plt.imshow(SampleImage)
    plt.colorbar()
    plt.show()
    
    print('\nDetected Reference Image')
    plt.figure()
    plt.imshow(ReferenceImage)
    plt.colorbar()
    plt.show()
       
    print('\nDetected Propagation Image')
    plt.figure()
    plt.imshow(PropagImage)
    plt.colorbar()
    plt.show()
    
    print('\nSubImage')
    plt.figure()
    plt.imshow(SubImage)
    plt.colorbar()
    plt.show()
    
    print('\ndxOF')
    plt.figure()
    plt.imshow(dxOF)
    plt.colorbar()
    plt.show()          
    print('\ndyOF')
    plt.figure()
    plt.imshow(dyOF)
    plt.colorbar()
    plt.show()          
    print('\nphiLarkin')
    plt.figure()
    plt.imshow(phiLarkin)
    plt.colorbar()
    plt.show()