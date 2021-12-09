'''
Authors
  - C. Selmi:  written in 2021
'''
import os
import numpy as np
import copy
from astropy.io import fits as pyfits
from LibSpecchi.type.modalAmplitude import ModalAmplitude
from LibSpecchi.type.modalBase import ModalBase
from LibSpecchi.type.commandHistory import CmdHistory
from Lib.tracking_number_folder import TtFolder

class IFMaker():
    '''
    '''

    def __init__(self, deformable_mirror, interferometer):
        """The constructor """
        self._dm = deformable_mirror
        self._interf = interferometer
        if deformable_mirror is not None:
            self._nActs = deformable_mirror.getNActs() #qualcosa del genere

        #acquisizione
        self._nPushPull = None
        self._template = None
        self._amplitudeTag = None
        self._amplitude = None
        self._cmdMatrixTag = None
        self._cmdMatrix = None
        self._actsVector = None
        self._tt_cmdH = None
        self._indexingList = None
        self._tt = None

        #analisi
        self._cube = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return '/Users/rm/Desktop/Arcetri/M4/Data/M4Data/OPTData/IFFunctions'

    def acquisition(self, n_push_pull, cmd_matrix_fits_file_name,
                        amplitude_fits_file_name,
                        shuffle=False, template=None):
        '''
        Performs the process of acquiring interferograms

        Parameters
        ----------
             n_push_pull: int
                          number of push pull for each mode
             amplitude_fits_file_name: string
                                       fits file name
                                       (Example = 'amp.fits' a vector of shape nActs)
             cmd_matrix_fits_file_name: string
                                        fits file name
                                        (Example = 'modalBase.fits' matrix of shape nActs x nActs)
        Other Parameters
        ----------------
             shuffle: boolean, optional
                      if not indicated, the function create the tidy command
                      history matrix
                      True for shuffle acquisition
             template: numpy array, optional
                       vector composed by 1 and -1
                       if not indicated, the function use the vector [1, -1, 1]

        Returns
        -------
                tt: string
                    tracking number of measurements made
        '''
        amplitude, cmd_matrix = self._readTypeFromFitsNameTag(amplitude_fits_file_name,
                                                              cmd_matrix_fits_file_name)

        self._nPushPull = n_push_pull
        if template is None:
            self._template = np.array([1, -1, 1])
        else:
            self._template = template
        self._amplitudeTag = amplitude_fits_file_name
        self._amplitude = amplitude
        self._cmdMatrixTag = cmd_matrix_fits_file_name
        self._cmdMatrix = cmd_matrix

        self._actsVector = np.arange(self._nActs)
        indexing_input = copy.copy(self._actsVector)
        dove, tt = TtFolder(self._storageFolder())._createFolderToStoreMeasurements()
        self._tt = tt

        cmdH = CmdHistory(self._nActs)
        if shuffle is False:
            command_history_matrix_to_apply, self._tt_cmdH = \
                    cmdH.tidyCommandHistoryMaker(self._actsVector,
                                                 amplitude,
                                                 cmd_matrix,
                                                 n_push_pull,
                                                 template)
        if shuffle is True:
            command_history_matrix_to_apply, self._tt_cmdH = \
                    cmdH.shuffleCommandHistoryMaker(self._actsVector,
                                                    amplitude,
                                                    cmd_matrix,
                                                    n_push_pull,
                                                    template)
        self._indexingList = cmdH.getIndexingList()
        self._saveInfo(dove)

        n_images = 1
        for i in range(command_history_matrix_to_apply.shape[1]):
            self._dm.setActsCommand(command_history_matrix_to_apply[:, i]) #vecchia setShape
            masked_image = self._interf.acquire_phasemap(n_images)
            file_name = 'image_%04d.fits' %i
            self._interf.save_phasemap(dove, file_name, masked_image)

        return tt

    def _saveInfo(self, folder):
        fits_file_name = os.path.join(folder, 'info.fits')
        header = pyfits.Header()
        header['NPUSHPUL'] = self._nPushPull
        header['TT_CMDH'] = self._tt_cmdH
        header['CMDMAT'] = self._cmdMatrixTag
        header['AMP'] = self._amplitudeTag
        header['NACTS'] = self._nActs
        pyfits.writeto(fits_file_name, self._indexingList, header)
        pyfits.append(fits_file_name, self._template)
        pyfits.append(fits_file_name, self._actsVector)


    def _readTypeFromFitsNameTag(self, amplitude_fits_file_name,
                                 cmd_matrix_fits_file_name):
        '''
        Parameters
        ----------
            amplitude_fits_file_name: string
                                     vector with mode amplitude fits file name
            cmd_matrix_fits_file_name: string
                                    matrix of mode commands fits file name

        Returns
        -------
            amplitude: numpy array
                    vector with mode amplitude
            cmd_matrix: numpy array [nActs x nModes]
                        matrix of mode commands
                        diagonal matrix in case of zonal commands
        '''
        ma = ModalAmplitude.loadFromFits(amplitude_fits_file_name)
        amplitude = ma.getModalAmplitude()

        mb = ModalBase.loadFromFits(cmd_matrix_fits_file_name)
        cmd_matrix = mb.getModalBase()
        return amplitude, cmd_matrix

    def analysis(self):
        cube = self._createCube()
        self._saveCube('Cube.fits')
        return cube

    def _createCube(self):
        '''
        Parameters
        ----------
                tt: string

        Returns
        -------
                cube = masked array [pixels, pixels, number of images]
                        cube from analysis
        '''

        cube_all_act = None
        where = self._indexReorganization(self._indexingList, self._actsVector,
                                          self._nPushPull)
        ampl_reorg = self._amplitudeReorganization(self._actsVector,
                                                   self._indexingList,
                                                   self._amplitude,
                                                   self._nPushPull)
        for i in range(self._actsVector.shape[0]):
            print(i)
            for k in range(self._nPushPull):
                p = self._nPushPull * i + k
                n = where[p]
                mis_amp = k* self._indexingList.shape[1] + n
                mis = k * self._indexingList.shape[1] * self._template.shape[0] \
                        + n * self._template.shape[0]

                name = 'image_%04d.fits' %mis
                file_name = os.path.join(self._storageFolder(), self._tt, name)
                image0 = self._interf.readImage4D(file_name) #creare questa funzione nell'interf prendendola da ic
                image_list = [image0]
                for l in range(1, self._template.shape[0]):
                    name = 'image_%04d.fits' %(mis+l)
                    file_name = os.path.join(self._storageFolder(), self._tt, name)
                    ima = self._interf.readImage4D(file_name)
                    image_list.append(ima)

                image = np.zeros((image0.shape[0], image0.shape[1]))
                for p in range(1, len(image_list)):
                    opd2add = image_list[p] * self._template[p] + image_list[p-1] * self._template[p-1]
                    master_mask2add = np.ma.mask_or(image_list[p].mask, image_list[p-1].mask)
                    if p==1:
                        master_mask = master_mask2add
                    else:
                        master_mask = np.ma.mask_or(master_mask, master_mask2add)
                    image += opd2add
                image = np.ma.masked_array(image, mask=master_mask)
                    #image = image + ima * self._template[l] #sbagliato
                img_if = image / (2 * ampl_reorg[mis_amp] * (self._template.shape[0] - 1))

                if_push_pull_kth = img_if

                if k == 0:
                    all_push_pull_act_jth = if_push_pull_kth
                else:
                    all_push_pull_act_jth = np.ma.dstack((all_push_pull_act_jth,
                                                          if_push_pull_kth))

            if self._nPushPull == 1:
                if_act_jth = all_push_pull_act_jth
            else:
                if_act_jth = np.ma.mean(all_push_pull_act_jth, axis=2)

            if cube_all_act is None:
                cube_all_act = if_act_jth
            else:
                cube_all_act = np.ma.dstack((cube_all_act, if_act_jth))
        self._cube = cube_all_act

        return self._cube

    def _indexReorganization(self, indexingList, actsVector, nPushPull):
        """ Returns the index position """
        indv = np.array(indexingList)
        where = []
        for ind in actsVector:
            for j in range(nPushPull):
                a = np.where(indv[j] == ind)
                where.append(a[0][0])
        return where

    def _amplitudeReorganization(self, indexing_input, indexing_list,
                                 amplitude, n_push_pull):
        '''
            Args:
                indexing_input = vector of selected modes for carrying out
                                 the influence functions

                indexing_list = tuple indicating how the modes were applied

                amplitude = amplitude of applied modes

            Returns:
                vect = vector (amp.shape x n_push_pull.shape) with the
                        amplitudes ordered in the same way as indexing_list

        '''
        where = []
        for i in indexing_input:
            for j in range(n_push_pull):
                a = np.where(indexing_list[j] == i)
                where.append(a)
        where = np.array(where)
        vect = np.zeros(amplitude.shape[0]*n_push_pull)

        for i in range(amplitude.shape[0]):
            for k in range(n_push_pull):
                p = n_push_pull * i + k
                indvect = where[p]+ indexing_input.shape[0] * k
                vect[indvect] = amplitude[i]
        return vect

    def _saveCube(self, cube_name):
        """
        Parameters
        ----------
                cube_name: string
                            name to save the cube
                            example 'Cube.fits'
        """
        dove = os.path.join(self._storageFolder(), self._tt)
        file_name = os.path.join(dove, cube_name)
        header = pyfits.Header()
        header['NPUSHPUL'] = self._nPushPull
        header['NACTS'] = self._nActs
        pyfits.writeto(file_name, self._cube.data, header)
        pyfits.append(file_name, self._cube.mask.astype(int))
        pyfits.append(file_name, self._amplitude)
        pyfits.append(file_name, self._actsVector)

    def getCube(self):
        '''
        Returns
        -------
                cube: masked array [pixels, pixels, number of images]
                    cube from analysis
        '''
        return self._cube

    @staticmethod
    def loadIFMaker(tt):
        """ Creates the object using information contained in Cube

        Parameters
        ----------
                fits_file_name: string
                                cube file name path

        Returns
        -------
                theObject: object
                            analyzerIFF class object
        """
        theObject = IFMaker(None, None)
        file_name = os.path.join(theObject._storageFolder(), tt, 'Cube.fits')
        header = pyfits.getheader(file_name)
        hduList = pyfits.open(file_name)
        cube = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        cmd_amplitude = hduList[2].data
        acts_vector = hduList[3].data
        try:
            n_push_pull = header['NPUSHPUL']
        except KeyError:
            n_push_pull = 1

        theObject._nActs = header['NACTS']
        theObject._actsVector = acts_vector
        theObject._amplitude = cmd_amplitude
        theObject._nPushPull = n_push_pull
        theObject._cube = cube
        return theObject
