'''
'''

import os
import configparser

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), 'config.conf'))

ROOT_FOLDER = config['path']['root_folder']

MODALAMPLITUDE_ROOT_FOLDER = os.path.join(ROOT_FOLDER, 'ModalAmplitude')
MODALBASE_ROOT_FOLDER = os.path.join(ROOT_FOLDER, 'ModalBase')
COMMANDHISTORY_ROOT_FOLDER = os.path.join(ROOT_FOLDER, 'CommandHistory')
IFFUNCTIONS_ROOT_FOLDER = os.path.join(ROOT_FOLDER, 'IFFunctions')
FLAT_ROOT_FOLD = os.path.join(ROOT_FOLDER, 'Flattening')

OPD_SERIES_ROOT_FOLD = os.path.join(ROOT_FOLDER, 'OPD_series')
LINEARITY_ROOT_FOLD = os.path.join(ROOT_FOLDER, 'Linearity')
