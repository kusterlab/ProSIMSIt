import glob
import json
import os
import shutil
from pathlib import Path

from oktoberfest import runner
from oktoberfest.data import Spectra
from oktoberfest.utils import Config, JobPool

from prosimsit.constants import PROSIT_CONFIG


def generate_oktoberfest_config(config, mzml_folder: Path, config_path: Path):
    """
    Generate a config.json file to be utilized for both Oktoberfest runs
    :param config: Dictionary of all config parameters generated from config.toml
    :param mzml_folder: Directory containing mzML files
    :param config_path: Path to config.json output file
    :return: None
    """
    oktoberfest_config = PROSIT_CONFIG.copy()
    oktoberfest_config['tag'] = config['general']['tmt_type']
    oktoberfest_config['inputs']['spectra'] = str(mzml_folder)
    oktoberfest_config['inputs']['search_results'] = config['inputs']['maxquant_results']
    oktoberfest_config['output'] = str(Path(config['general']['output']) / 'oktoberfest_1_out')
    oktoberfest_config['models'] = {'intensity': config['prosit']['intensity_model'],
                                    'irt': config['prosit']['irt_model']}
    oktoberfest_config['prediction_server'] = config['prosit']['prediction_server']
    oktoberfest_config['numThreads'] = config['general']['threads']
    oktoberfest_config['thermoExe'] = None
    if config['prosit']['ssl']:
        oktoberfest_config['ssl'] = True
    if config['prosit']['ms_analyzer'] == 'ot':
        oktoberfest_config['massTolerance'] = 20
        oktoberfest_config['unitMassTolerance'] = 'ppm'
    elif config['prosit']['ms_analyzer'] == 'it':
        oktoberfest_config['massTolerance'] = 0.35
        oktoberfest_config['unitMassTolerance'] = 'da'
    elif config['prosit']['ms_analyzer'] == 'manual':
        oktoberfest_config['massTolerance'] = config['prosit']['mass_tolerance']
        if config['prosit']['tolerance_unit'] not in ['da', 'ppm']:
            raise ValueError(
                f"Unknown tolerance unit: {config['prosit']['tolerance_unit']}. "
                f"Use 'da' for dalton when rescoring ITMS data,"
                f"or 'ppm' for parts per million when rescoring FTMS data."
            )
        oktoberfest_config['unitMassTolerance'] = config['prosit']['tolerance_unit']
    else:
        raise ValueError(
            f"Unknown mass analyzer: {config['prosit']['ms_analyzer']}. Use 'it' for ion trap or 'ot' for orbitrap.")

    with open(config_path, 'w') as outfile:
        json.dump(oktoberfest_config, outfile, indent=4, )


def prepare_second_oktoberfest_run(mzml_dir, oktoberfest_config_path, msms_dir, output_dir):
    """
    This function prepares the second Oktoberfest run by copying the results from the first run to the second run.
    :param mzml_dir: Folder containing mzML files
    :param oktoberfest_config_path: Path to the config.json file generated from generate_oktoberfest_config()
    :param msms_dir: Folder containing MaxQuant msms.txt file
    :param output_dir: Path to the output directory
    :return: Oktoberfest Config object containing all oktoberfest-related configurations
    """
    conf = Config()
    conf.read(oktoberfest_config_path)
    conf.check()
    original_output_dir = conf.data['output']
    conf.inputs['search_results'] = msms_dir
    conf.data['output'] = output_dir / 'oktoberfest_2_out'
    conf.inputs['spectra'] = mzml_dir
    conf.inputs['spectra_type'] = 'mzml'

    def copy_files_with_pattern(source_dir, dest_dir, pattern):
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        search_pattern = os.path.join(source_dir, pattern)

        files = glob.glob(search_pattern)

        for file_path in files:
            shutil.copy(file_path, dest_dir)
            print(f"Copied {file_path} to {dest_dir}")

    source_directory = original_output_dir + '/results/'
    destination_directory = conf.output / 'results/'
    file_pattern = '*.txt'

    copy_files_with_pattern(source_directory, destination_directory, file_pattern)

    # copy progress files as well to indicate that CE calib is already done
    source_directory = original_output_dir + '/proc/'
    destination_directory = conf.output / 'proc/'
    file_pattern = 'ce_calib*'
    copy_files_with_pattern(source_directory, destination_directory, file_pattern)
    return conf


def preprocess_spectra_files(conf):
    """
    Wrapper to apply oktoberfest preprocessing steps to spectra files.
    :param conf: Config object for Oktoberfest generated from prepare_second_oktoberfest_run()
    :return: List of preprocessed spectra files
    """
    spectra_files = [Path(f) for f in glob.iglob(str(conf.inputs['spectra'] / '*'))]
    spectra_files = runner._preprocess(spectra_files, conf)
    return spectra_files


def annotate_library(spectra_files, conf):
    """
    Wrapper to apply oktoberfest library annotation steps.
    :param spectra_files: List of preprocessed spectra files
    :param conf: Config object for Oktoberfest generated from prepare_second_oktoberfest_run()
    :return: None
    """
    if conf.num_threads > 1:
        processing_pool = JobPool(processes=conf.num_threads)
        for spectra_file in spectra_files:
            processing_pool.apply_async(runner._annotate_and_get_library, [spectra_file, conf])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            runner._annotate_and_get_library(spectra_file, conf)


def generate_pred_files(conf):
    """
    Generate predictions for all hdf5 files prepared for Prosit.
    :param conf: Config object for Oktoberfest generated from prepare_second_oktoberfest_run()
    :return: None
    """
    spectra_files_str = [f for f in glob.iglob(str(conf.output / 'data' / '*'))]
    for f in spectra_files_str:
        library = Spectra.from_hdf5(f)
        with open(conf.output / 'results/' / (f.split('/')[-1].split('.')[0] + '_ce.txt'), 'r') as file:
            content = file.read()
            best_ce = int(content)
        library.obs['COLLISION_ENERGY'] = best_ce
        library.write_as_hdf5(f.replace('.hdf5', '.pred.hdf5'))


def calculate_featuers(spectra_files, conf):
    """
    Function to calculate Percolator features from Prosit predictions.
    :param spectra_files: List of preprocessed spectra files
    :param conf: Config object for Oktoberfest generated from prepare_second_oktoberfest_run()
    :return: None
    """
    if conf.num_threads > 1:
        processing_pool = JobPool(processes=conf.num_threads)
        for spectra_file in spectra_files:
            processing_pool.apply_async(runner._calculate_features, [spectra_file, conf])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            runner._calculate_features(spectra_file, conf)
