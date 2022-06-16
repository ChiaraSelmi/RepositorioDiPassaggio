'''
Authors
  - C. Selmi:  written in 2019 and 2022
'''
import numpy as np
import os
import glob
from astropy.io import fits
from LibSpecchi.configuration import config

class Analyser():
    
    def __init__(self):
        pass
    
    def linearity(self, tn):
        fold_for_meas = os.path.join(config.LINEARITY_ROOT_FOLD, tn)
        tt_list = os.listdir(fold_for_meas)
        del(tt_list[0])
        cube = None
        for tt in tt_list:
            path = glob.glob(os.path.join(fold_for_meas, tt, 'Cube.fits'))
            hduList = fits.open(path)
            cube_amp = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
            if cube is None:
                cube = cube_amp
            else:
                cube = np.ma.dstack((cube, cube_amp))
        return cube
