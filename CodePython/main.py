#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:37:59 2020

@author: quenot
"""
import datetime
import time
import os
from InputOutput.pagailleIO import saveEdf
from Experiment_multiprocessing import Experiment


def main(experiment_name, sample_sampling, nb_exp_points, margin,
         simulation_type, filepath, save=True, multiprocessing=True,
         cpus=None):
    """Main of the simulation code

    Parameters
    ----------
    experiment_name : str
        The name of the experiment in the .xml file
    sample_sampling : int
        Oversampling factor of the sample
    nb_exp_points : int
        Number of different membrane positions
    margin : int
        Prevent aliasing in Fresnel by extending images
    simulation_type : str
        The propagation method, must be Fresnel or RayT
    filepath : str
        The filepath to save the simulated images
    save : bool, optional
        Whether to save the output of the experiment in filepath, by default
        True
    multiprocessing : bool, optional
        Whether to process multiple energies in parallel, by default True
    cpus : int or None, optional
        The number of cores to use in the multiprocessing, by default None =
        All

    Raises
    ------
    ValueError
        If the given filepath is not found
    ValueError
        If simulation_type is neither Fresnel nor RayT
    """
    time0 = time.time()  # timer for computation
    exp_dict = {}

    # Define experiment
    exp_dict['experimentName'] = experiment_name
    exp_dict['filepath'] = filepath
    exp_dict['sampleSampling'] = sample_sampling
    exp_dict['nbExpPoints'] = nb_exp_points
    exp_dict['margin'] = margin
    exp_dict['simulation_type'] = simulation_type
    exp_dict['Multiprocessing'] = multiprocessing
    exp_dict['CPUs'] = cpus

    if save and not os.path.isdir(filepath):
        raise ValueError('Path not found')

    # ************************************************************************
    now = datetime.datetime.now()
    exp_dict['expID'] = now.strftime("%Y%m%d-%H%M%S")  # define experiment ID

    experiment = Experiment(exp_dict)

    for point_num in range(exp_dict['nbExpPoints']):
        experiment.myMembrane.myGeometry = []
        experiment.myMembrane. \
            getMyGeometry(experiment.studyDimensions,
                          experiment.myMembrane.membranePixelSize,
                          experiment.sampling, point_num,
                          exp_dict['nbExpPoints'])
        print("\n\nINITIALIZING EXPERIMENT PARAMETERS AND GEOMETRIES")
        print("\n\n*************************")
        print("Calculations point", point_num)
        if exp_dict['simulation_type'] == "Fresnel":
            sample_image_tmp, reference_image_tmp, propag_image_tmp, white = \
                experiment.computeSampleAndReferenceImages(point_num)
        elif exp_dict['simulation_type'] == "RayT":
            sample_image_tmp, reference_image_tmp, propag_image_tmp, white = \
                experiment.computeSampleAndReferenceImagesRT(point_num)
        else:
            raise ValueError(f"Simulation Type ({exp_dict['simulation_type']})\
            must be either Fresnel or RayT")
        nbin = len(sample_image_tmp)

        if point_num == 0:
            exp_path_en = []
            if exp_dict['simulation_type'] == "Fresnel":
                exp_images_filepath = exp_dict['filepath'] + \
                    'Fresnel_'+str(exp_dict['expID'])+'/'
            if exp_dict['simulation_type'] == "RayT":
                exp_images_filepath = exp_dict['filepath'] + \
                    'RayTracing_'+str(exp_dict['expID'])+'/'
            if save:
                os.mkdir(exp_images_filepath)
                os.mkdir(exp_images_filepath+'membraneThickness/')
            thresholds = experiment.myDetector.myBinsThresholds.copy()
            thresholds.insert(0, experiment.mySource.mySpectrum[0][0])
            for ibin in range(nbin):
                binstart = '%2.2d' % thresholds[ibin]
                binend = '%2.2d' % thresholds[ibin+1]
                exp_path_en.append(f'{exp_images_filepath}{binstart}_\
                    {binend}kev/')
                if len(thresholds)-1 == 1:
                    exp_path_en = [exp_images_filepath]
                if save:
                    os.mkdir(exp_path_en[ibin])
                    os.mkdir(exp_path_en[ibin]+'ref/')
                    os.mkdir(exp_path_en[ibin]+'sample/')
                    os.mkdir(exp_path_en[ibin]+'propag/')

        txt_point = '%2.2d' % point_num
        if save:
            saveEdf(experiment.myMembrane.myGeometry[0], exp_images_filepath +
                    'membraneThickness/'+exp_dict['experimentName'] +
                    '_sampling'+str(exp_dict['sampleSampling'])+'_' +
                    str(point_num)+'.edf')
            for ibin in range(nbin):
                saveEdf(sample_image_tmp[ibin], exp_path_en[ibin] +
                        'sample/sampleImage_'+str(exp_dict['expID'])+'_' +
                        txt_point+'.edf')
                saveEdf(reference_image_tmp[ibin], exp_path_en[ibin] +
                        'ref/ReferenceImage_'+str(exp_dict['expID'])+'_' +
                        txt_point+'.edf')
                if point_num == 0:
                    saveEdf(propag_image_tmp[ibin], exp_path_en[ibin] +
                            'propag/PropagImage_'+str(exp_dict['expID'])+'_' +
                            '.edf')
                    saveEdf(white[ibin], exp_path_en[ibin]+'white_' +
                            str(exp_dict['expID'])+'_'+'.edf')
    if save:
        experiment.saveAllParameters(time0, exp_dict)

    print("\nfini")


if __name__ == '__main__':
    main("SIMAP_SpheresInTube", 2, 1, 10, 'RayT', '/Users/chris/Documents/')
