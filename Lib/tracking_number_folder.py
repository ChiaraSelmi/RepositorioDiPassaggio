'''
@author: cs
'''

import os
from Lib.timestamp import Timestamp

class TtFolder():
    """ Create a new folder using the generation of the tracking number """

    def __init__(self, store_in_folder):
        """The constructor """
        self._rootStoreFolder = store_in_folder


    def _createFolderToStoreMeasurements(self):
        """
        returns:
            dove = all path generated
            tt = only tracking number
        """
        tt = Timestamp.now()
        dove = os.path.join(self._rootStoreFolder, tt)
        if os.path.exists(dove):
            raise OSError('Directory %s exists' % dove)
        else:
            os.makedirs(dove)
        return dove, tt
