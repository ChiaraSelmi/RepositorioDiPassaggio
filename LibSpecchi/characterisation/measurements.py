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
from LibSpecchi.convertWFToDmCommand import Converter

class MeasurementAcquisition():

    def __init__(self, deformable_mirror, interferometer, tn_iff):
        self.dm = deformable_mirror
        self.interf = interferometer
        self.converter = Converter(tn_iff)

    def flattening(self, cmd=None):
        if cmd is None:
            wf = self.interf.wavefront()
            command = self.converter.fromWfToDmCommand(wf)
        else:
            command = cmd

        pos = self.dm.get_shape()
        cmd_to_apply = pos + command
        maxComm = max(abs(cmd_to_apply))
        print('max command= %f' % maxComm)

        if maxComm>config.MAX_COMMAND_TO_APPLY:
            raise OSError('Actuator command too large')
        else:
            self.dm.set_shape(-cmd_to_apply)

        self._commandToApply = cmd_to_apply

    def closeLoop(self, n_meas):
        tt_list = []
        for i in range(0, n_meas):
            fold_for_meas = config.FLAT_ROOT_FOLD
            dove, tt = TtFolder(fold_for_meas).createFolderToStoreMeasurements()
            wf = self.interf.wavefront()
            cmd = self.converter.fromWfToDmCommand(wf)
            fits.writeto(dove + 'imgstart.fits', wf.data)
            fits.append(dove + 'imgstart.fits', wf.mask.astype(int))
            fits.writeto(dove + 'flatDeltaCommand.fits', cmd)

            self.flattening(cmd)
            wf = self.interf.wavefront()
            fits.writeto(dove + 'imgflat.fits', wf.data)
            fits.append(dove + 'imgflat.fits', wf.mask.astype(int))

            fits.writeto(dove + 'flatCommand.fits', self._commandToApply)
            tt_list.append(tt)
        return tt_list

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
        Nella cartella finale sono comprese una cartella per ogni valore di amp
        presente nel vettore in ingresso dato che per ogni amp si usa l'iffMaker

        Parameters
        ----------
        cmd_matrix_tag: string
            tag for command matrix containing 1 on the diagonal for the actuators
            to be used
        vector_of_amplitude_for_meas: numpy array
            vector containing the amplitude values to be given to the actuators
            to make the measurement

        Returns
        -------
        tt = string
            tracking number of measurements
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
        for amp in amp_list_tag:
            tn = iff.acquisitionAndAnalysis(cmd_matrix_tag,
                                            amp, shuffle=False,
                                            template=None)
            source = os.path.join(iff._storageFolder(), tn)
            destination = dove
            shutil.move(source, destination)
        fits.writeto(os.path.join(dove, 'amplitude.fits'),
                     vector_of_amplitude_for_meas, overwrite=True)
        return tt

