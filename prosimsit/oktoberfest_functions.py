import glob
import os
import shutil
from pathlib import Path

from oktoberfest import runner
from oktoberfest.data import Spectra
from oktoberfest.utils import Config, JobPool


def prepare_second_oktoberfest_run(mzml_dir, oktoberfest_config_path, dir_msms, output_dir):
    conf = Config()
    conf.read(oktoberfest_config_path)
    conf.check()
    original_output_dir = conf.data['output']
    conf.inputs['search_results'] = dir_msms
    conf.data['output'] = output_dir / 'oktoberfest_2_out'
    conf.inputs['spectra'] = mzml_dir
    conf.inputs['spectra_type'] = 'mzml'

    def copy_files_with_pattern(source_dir, dest_dir, pattern):
        # Create destination directory if it doesn't exist
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        # Create a pattern to search for files
        search_pattern = os.path.join(source_dir, pattern)

        # Find files that match the pattern
        files = glob.glob(search_pattern)

        # Copy each file to the destination directory
        for file_path in files:
            shutil.copy(file_path, dest_dir)
            print(f"Copied {file_path} to {dest_dir}")

    source_directory = original_output_dir + '/results/'
    destination_directory = conf.output / 'results/'
    file_pattern = '*.txt'  # Pattern for all .txt files

    copy_files_with_pattern(source_directory, destination_directory, file_pattern)

    # copy progress files as well to indicate that CE calib is already done
    source_directory = original_output_dir + '/proc/'
    destination_directory = conf.output / 'proc/'
    file_pattern = 'ce_calib*'  # Pattern for all .txt files
    copy_files_with_pattern(source_directory, destination_directory, file_pattern)
    return conf


def preprocess_spectra_files(conf):
    spectra_files = [Path(f) for f in glob.iglob(str(conf.inputs['spectra'] / '*'))]
    spectra_files = runner._preprocess(spectra_files, conf)
    return spectra_files


def annotate_library(spectra_files, conf):
    if conf.num_threads > 1:
        processing_pool = JobPool(processes=conf.num_threads)
        for spectra_file in spectra_files:
            processing_pool.apply_async(runner._annotate_and_get_library, [spectra_file, conf])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            runner._annotate_and_get_library(spectra_file, conf)


def generate_pred_files(conf):
    spectra_files_str = [f for f in glob.iglob(str(conf.output / 'data' / '*'))]
    for f in spectra_files_str:
        library = Spectra.from_hdf5(f)
        with open(conf.output / 'results/' / (f.split('/')[-1].split('.')[0] + '_ce.txt'), 'r') as file:
            content = file.read()
            best_ce = int(content)
        library.obs['COLLISION_ENERGY'] = best_ce
        library.write_as_hdf5(f.replace('.hdf5', '.pred.hdf5'))


def calculate_featuers(spectra_files, conf):
    if conf.num_threads > 1:
        processing_pool = JobPool(processes=conf.num_threads)
        for spectra_file in spectra_files:
            processing_pool.apply_async(runner._calculate_features, [spectra_file, conf])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            runner._calculate_features(spectra_file, conf)