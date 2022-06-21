'''
Authors
  - C. Selmi:  written in 2019 and 2022
'''
import numpy as np
import os
#import glob
from astropy.io import fits
from LibSpecchi.configuration import config
from LibSpecchi.influenceFunctionsMaker import IFMaker

class Analyser():

    def __init__(self):
        pass

    def linearity(self, tn):
        fold_for_meas = os.path.join(config.LINEARITY_ROOT_FOLD, tn)
        tt_list = os.listdir(fold_for_meas)
        #del(tt_list[0])
        del(tt_list[-1])
        cube = None
        for tt in tt_list:
            iff = IFMaker.loadAnalyzerFromIFMaker(os.path.join(tn, tt),
                                                  config.LINEARITY_ROOT_FOLD)
            cube_amp = iff.getCube()
            if cube is None:
                cube = cube_amp
            else:
                cube = np.ma.dstack((cube, cube_amp))

        rms_list=[]
        for i in range(cube.shape[2]):
            rms = np.ma.std(cube[:,:,i])
            rms_list.append(rms)
        rms_vect = np.array(rms_list)

        hduList = fits.open(os.path.join(fold_for_meas, 'amplitude.fits'))
        amp = hduList[0].data

        linearityMatrix = np.reshape(rms_vect,
                                     (amp.size, np.int(cube.shape[2]/amp.size)))
        return linearityMatrix
