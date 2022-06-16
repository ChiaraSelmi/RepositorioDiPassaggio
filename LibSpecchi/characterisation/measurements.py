'''
Authors
  - C. Selmi:  written in 2019 and 2022
'''
import os
import numpy as np
import time
from astropy.io import fits
import shutil
from LibSpecchi.ground import zernike
from LibSpecchi.ground.timestamp import Timestamp
from LibSpecchi.influenceFunctionsMaker import IFMaker
from LibSpecchi.ground.tracking_number_folder import TtFolder
from LibSpecchi.configuration import config
from LibSpecchi.type.modalBase import ModalBase
from LibSpecchi.type.modalAmplitude import ModalAmplitude

class MeasurementAcquisition():

    def __init__(self, deformable_mirror, interferometer, tn_iff):
        self.dm = deformable_mirror
        self.interf = interferometer
        self.an = IFMaker.loadAnalyzerFromIFMaker(tn_iff)

    def opticalMonitoring(self, n_images, delay):
        '''
        Parameters
        ----------
        n_images: int
            number of images to acquire
        delay: int [s]
            waiting time (in seconds) between two image acquisitions

        Returns
        ------
        tt: string
            tracking number of measurements
        '''
        fold_for_meas = config.OPD_SERIES_ROOT_FOLD
        dove, tt = TtFolder(fold_for_meas).createFolderToStoreMeasurements()

        zer_list = []
        t0 = time.time()
        for i in range(n_images):
            ti = time.time()
            dt = ti - t0
            masked_ima = self.interf.wavefront(1)
            name = Timestamp.now() + '.fits'
            fits_file_name = os.path.join(dove, name)
            fits.writeto(fits_file_name, masked_ima.data)
            fits.append(fits_file_name, masked_ima.mask.astype(int))

            coef, mat = zernike.zernikeFit(masked_ima, np.arange(10) + 1)
            vect = np.append(dt, coef)
            zer_list.append(vect)

            fits_file_name = os.path.join(dove, 'zernike.fits')
            fits.writeto(fits_file_name, np.array(zer_list), overwrite=True)

            time.sleep(delay)

    def linearity(self, cmd_matrix_tag, vector_of_amplitude_for_meas):
        '''
        la cmd matrix se la costruisce l'utente mettendo 1 agli attuatori che vuole
        misurare. L'ampiezza la creo io.
        '''
        mb = ModalBase.loadFromFits(cmd_matrix_tag)
        cmd_matrix = mb.getModalBase()
        ma = ModalAmplitude()
        amp_list_tag = []
        for i in range(vector_of_amplitude_for_meas.size):
            amp = vector_of_amplitude_for_meas[i]
            single_amp = np.zeros(cmd_matrix.shape[0]) + amp
            tag = 'amp%.02f' %amp
            ma.saveAsFits(tag, single_amp)
            amp_list_tag.append(tag)


        fold_for_meas = config.LINEARITY_ROOT_FOLD
        dove, tt = TtFolder(fold_for_meas).createFolderToStoreMeasurements()
        iff = IFMaker(self.interf, self.dm)
        n_rep = 1
        for amp in amp_list_tag:
            tn = iff.acquisitionAndAnalysis(n_rep, cmd_matrix_tag,
                                            amp, shuffle=False,
                                            template=None)
            source = os.path.join(iff._storageFolder(), tn)
            destination = dove
            shutil.move(source, destination)
        return tt
