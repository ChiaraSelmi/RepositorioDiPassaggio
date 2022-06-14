'''
Authors
  - C. Selmi: written in 2019
'''
import os
import copy
import logging
from astropy.io import fits as pyfits
import h5py
import numpy as np
from LibSpecchi.ground.tracking_number_folder import TtFolder
from LibSpecchi.configuration import config


class CmdHistory():
    '''
    Class to create and manage the command history matrix

    HOW TO USE IT::

        from m4.type.commandHistory import CmdHistory
        cmdH = CmdHistory()
    '''

    def __init__(self, nActs):
        """The constructor """
        self._logger = logging.getLogger('CMD_HISTORY:')
        self._nActs = nActs
        self._modeVector = None
        self._template = None
        self._cmdMatrix = None
        self._ampVect = None
        self._nRepetitions = None
        #costruite dati gli input
        self._cmdSequence = None
        self._cmdHToApply = None

    def getCommandHistory(self):
        '''
        Returns
        -------
        ccmdHToApply: numpy array
                    command history matrix to apply
        '''
        return self._cmdHToApply

    def getIndexingList(self):
        '''
        Returns
        -------
            indexingList: list
                        list of modes/actuators used to create the command history matrix
        '''
        return self._indexingList

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save data"""
        # '/Users/rm/Desktop/Arcetri/M4/Data/M4Data/OPTData/CommandHistory'
        return config.COMMANDHISTORY_ROOT_FOLDER


    def tidyCommandHistoryMaker(self, mode_vector, amp_vector,
                                cmd_matrix, n_rep, template=None):
        '''
        Parameters
        ----------
            modesVector: numpy array
                         Mode or actuator index vector to apply
            ampVector: numpy array
                        amplitude mode vector
            cmdMatrix: numpy array [nActs x nModes]
                        mode command matrix
                        diagonal matrix in case of zonal commands
            n_rep: int
                        number of repetitions for each measurement mode 
        Other Parameters
        ----------
            template: numpy array  , optional
                    vector composed by 1 and -1
                    (es. np.array([1, -1, 1]))
        Returns
        -------
             matrixToApply: numpy array [nAct, nModes x n_rep x template.size]
                            tidy command history
             tt: string
                 tracking number
        '''
        self._ampVect = amp_vector
        self._cmdSequence = self._tidyCmdSequence(mode_vector, n_rep, cmd_matrix)
        aa = np.arange(self._cmdSequence.shape[1])
        bb = np.tile(amp_vector, n_rep)
        zipped = zip(aa, bb)
        self._cmdHToApply = self._cmdHistoryToApply(zipped, template)
        tt = self.saveInfo(0)
        self._logger.info('Creation of the ordered commandHistoryMatrix %s', tt)
        print(tt)

        return self._cmdHToApply, tt


    def shuffleCommandHistoryMaker(self, mode_vector, amp_vector,
                                   cmd_matrix, n_rep, template=None):
        '''
        Parameters
        ----------
            modesVector: numpy array
                         Mode or actuator index vector to apply
            ampVector: numpy array
                        amplitude mode vector
            cmdMatrix: numpy array [nActs x nModes]
                        mode command matrix
                        diagonal matrix in case of zonal commands
            n_rep: int
                        number of repetitions for each measurement mode 
        Other Parameters
        ----------
            template: numpy array  , optional
                    vector composed by 1 and -1
                    (es. np.array([1, -1, 1]))
        Returns
        -------
             matrixToApply: numpy array [nAct, nModes x n_rep x template.size]
                            shuffle command history
             tt: string
                 tracking number
        '''
        self._ampVect = amp_vector
        self._cmdSequence, self._indexingList = self._shuffleCmdSequence(
            mode_vector, n_rep, cmd_matrix)
        zipped = self._zippedAmplitude(amp_vector)
        self._cmdHToApply = self._cmdHistoryToApply(zipped, template)
        tt = self.saveInfo(0)
        self._logger.info('Creation of the shuffle commandHistoryMatrix %s', tt)
        print(tt)

        return self._cmdHToApply, tt

    def _shuffleCmdSequence(self, mode_vector, n_rep, cmd_matrix):
        ''' Usa il numero di ripetioni e crea la sequenza dei comandi
        '''
        self._modeVector = copy.copy(mode_vector)
        self._nRepetitions = n_rep
        self._cmdMatrix = cmd_matrix

        n_frame = mode_vector.size * n_rep
        matrix_for_seq = np.zeros((self._nActs, n_frame))

        indexingList = []
        for j in range(n_rep):
            np.random.shuffle(mode_vector)
            indexingList.append(list(mode_vector))

            cmdList = []
            for i in mode_vector:
                cmd = cmd_matrix[:, i]
                cmdList.append(cmd)

            for i in range(len(cmdList)):
                k = len(cmdList)*j + i
                matrix_for_seq.T[k] = cmdList[i]

        cmdSequence = matrix_for_seq
        indexingVec = np.array(indexingList)

        return cmdSequence, indexingVec


    def _tidyCmdSequence(self, mode_vector, n_rep, cmd_matrix):
        ''' Usa il numero di ripetioni e crea la sequenza dei comandi
        '''
        self._modeVector = copy.copy(mode_vector)
        self._nRepetitions = n_rep
        self._cmdMatrix = cmd_matrix
        indList = []
        for i in range(n_rep):
            indList.append(mode_vector)
        self._indexingList = np.array(indList)
        #self._indexingList= np.tile(mode_vector, n_push_pull)

        n_frame = mode_vector.size * n_rep
        matrix_for_seq = np.zeros((self._nActs, n_frame))

        cmdList = []
        for i in mode_vector:
            cmd = cmd_matrix[:, i]
            cmdList.append(cmd)

        for j in range(n_rep):
            for i in range(len(cmdList)):
                k = len(cmdList)*j + i
                matrix_for_seq.T[k] = cmdList[i]

        cmdSequence = matrix_for_seq

        return cmdSequence


    def _amplitudeReorganization(self, indexing_input, indexing_list,
                                 amplitude, n_push_pull):
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
                indvect = where[p][0][0]+ indexing_input.shape[0] * k
                vect[indvect] = amplitude[i]

        return vect


    def _zippedAmplitude(self, amplitude):
        aa = np.arange(self._cmdSequence.shape[1])
        reorganized_amplitude = \
                self._amplitudeReorganization(self._modeVector,
                                              self._indexingList,
                                              amplitude,
                                              self._nRepetitions)
        #zipped= np.dstack((aa, reorganized_amplitude))
        zipped = zip(aa, reorganized_amplitude)
        return zipped

    def _cmdHistoryToApply(self, zipped, template=None):
        ''' Usa l'informazione sul template per creare la sequenza di comandi
        da applicare allo specchio (il numero di ripetioni è già stato usato
        in precedenza per creare 
        '''
        matrix_with_amp = self._cmdSequence
        for i, amp in zipped:
            matrix_with_amp.T[i] = matrix_with_amp.T[i] * amp

        if template is None:
            template = np.array((1, -1, 1))
        else:
            template = template

        matrix_to_apply = np.zeros((self._nActs,
                                    self._cmdSequence.shape[1] *
                                    template.shape[0]))

        for i in range(self._cmdSequence.shape[1]):
            j = template.shape[0] * i
            for k in range(template.shape[0]):
                matrix_to_apply.T[j+k] = matrix_with_amp.T[i]* template[k]

        return matrix_to_apply


    def saveInfo(self, fits_or_h5):
        """ Save the data in fits format

        Returns
        -------
        tt: string
            tracking number
        """
        store_in_folder = CmdHistory._storageFolder()
        dove, tt = TtFolder(store_in_folder)._createFolderToStoreMeasurements()

        if fits_or_h5 == 0:
            fits_file_name = os.path.join(dove, 'info.fits')
            header = pyfits.Header()
            header['NREP'] = self._nRepetitions
            pyfits.writeto(fits_file_name, self._modeVector, header)
            pyfits.append(fits_file_name, self._indexingList, header)
            pyfits.append(fits_file_name, self._cmdMatrix, header)
            pyfits.append(fits_file_name, self._cmdHToApply, header)
            pyfits.append(fits_file_name, self._ampVect, header)
        else:
            fits_file_name = os.path.join(dove, 'info.h5')
            hf = h5py.File(fits_file_name, 'w')
            hf.create_dataset('dataset_1', data=self._modeVector)
            hf.create_dataset('dataset_2', data=self._indexingList)
            hf.create_dataset('dataset_3', data=self._cmdHToApply)
            hf.create_dataset('dataset_4', data=self._ampVect)
            hf.attrs['NREP'] = self._nRepetitions
            hf.close()
        return tt

    @staticmethod
    def load(tt, fits_or_h5=0):
        """ Creates the object from the info.fits file located in tt

        Parameters
        ----------
        tt: string
            tracking number

        Returns
        -------
        theObject: object
                command history class object
        """
        theObject = CmdHistory(None)
        theObject._tt = tt
        store_in_folder = CmdHistory._storageFolder()
        folder = os.path.join(store_in_folder, tt)
        if fits_or_h5 == 0:
            additional_info_fits_file_name = os.path.join(folder, 'info.fits')
            header = pyfits.getheader(additional_info_fits_file_name)
            hduList = pyfits.open(additional_info_fits_file_name)
            theObject._modeVector = hduList[0].data
            theObject._indexingList = hduList[1].data
            theObject._cmdMatrix = hduList[2].data
            theObject._cmdHToApply = hduList[3].data
            theObject._ampVect = hduList[4].data
            try:
                theObject._nRepetitions = header['NREP']
            except KeyError:
                theObject._nRepetitions = 1
        else:
            file_name = os.path.join(folder, 'info.h5')
            hf = h5py.File(file_name, 'r')
            hf.keys()
            data1 = hf.get('dataset_1')
            data2 = hf.get('dataset_2')
            data3 = hf.get('dataset_3')
            data4 = hf.get('dataset_4')
            theObject._nRepetitions = hf.attrs['NREP']
            theObject._modeVector = np.array(data1)
            theObject._indexingList = np.array(data2)
            theObject._cmdHToApply = np.array(data3)
            theObject._ampVect = np.array(data4)
            hf.close()
        return theObject
