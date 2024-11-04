import shutil
import logging

import pandas as pd

from simsi_transfer import maxquant as mq
from simsi_transfer import simsi_output
from simsi_transfer import evidence

logger = logging.getLogger(__package__ + "." + __file__)

def prepare_simsi_files(output_folder, maxquant_folder):
    """
    This function prepares the files for the SIMSI-Transfer workflow.
    """
    if (output_folder / 'simsi_input' / 'msmsScans.txt').is_file():
        logger.info('msmsScans.txt already exists, reusing it')
    else:
        shutil.copy(maxquant_folder / 'msmsScans.txt', output_folder / 'simsi_input' / 'msmsScans.txt')

    if (output_folder / 'simsi_input' / 'allPeptides.txt').is_file():
        logger.info('allPeptides.txt already exists, reusing it')
    else:
        shutil.copy(maxquant_folder / 'allPeptides.txt', output_folder / 'simsi_input' / 'allPeptides.txt')
    if (output_folder / 'simsi_input' / 'evidence.txt').is_file():
        logger.info('evidence.txt already exists, reusing it')
    else:
        shutil.copy(maxquant_folder / 'evidence.txt', output_folder / 'simsi_input' / 'evidence.txt')


def build_evidence(path_to_simsi_msms, mq_txt_folder, output_folder):
    if (output_folder / 'evidence.txt').is_file():
        logger.info('evidence.txt already exists, reusing it')
        return
    mq_txt_folders = [mq_txt_folder]
    msms_simsi = pd.read_csv(path_to_simsi_msms, sep='\t')
    logger.info(f'successfully read msms_simsi')

    evidence_mq = mq.read_evidence_txt(mq_txt_folder)
    allpeptides_mq = mq.read_allpeptides_txt(mq_txt_folder)
    plex = mq.get_plex(mq_txt_folders)

    logger.info(f'Starting SIMSI-Transfer evidence.txt building')
    evidence_simsi = evidence.build_evidence(msms_simsi, evidence_mq, allpeptides_mq, plex)
    evidence_simsi.to_csv(output_folder / 'evidence.txt', sep='\t', index=False, na_rep='NaN')