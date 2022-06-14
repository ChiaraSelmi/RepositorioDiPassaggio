'''
Authors
  - C. Selmi:  written in
'''

import os
from astropy.io import fits
import numpy as np
import h5py

def interf_save_phasemap(location, file_name, masked_image):
    """
    Parameters
    ----------
    location: string
        measurement file path
    file_name: string
        measuremnet fits file name
    masked_image: numpy masked array
        data to save
    """
    fits_file_name = os.path.join(location, file_name)
    fits.writeto(fits_file_name, masked_image.data)
    fits.append(fits_file_name, masked_image.mask.astype(int))

def interf_readImage4D4020(file_name):
    """
    Parameters
    ----------
        h5filename: string
             path of h5 file to convert
    Returns
    -------
            ima: numpy masked array
                 masked array image
    """
    file = h5py.File(file_name, 'r')
    genraw = file['measurement0']['genraw']['data']
    data = np.array(genraw)
    mask = np.zeros(data.shape, dtype=np.bool)
    mask[np.where(data == data.max())] = True
    ima = np.ma.masked_array(data * 632.8e-9, mask=mask)
    return ima

def interf_readImage(file_name):
    '''
    Parameters
    ----------
    file_name: string
        fits file path name of image to read
    '''
    hduList = fits.open(file_name)
    masked_ima = np.ma.masked_array(hduList[0].data,
                                    hduList[1].data.astype(bool))
    return masked_ima
