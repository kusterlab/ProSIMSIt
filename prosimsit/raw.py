from pathlib import Path
from typing import List, Optional

from simsi_transfer.thermo_raw import convert_raw_mzml_batch

from prosimsit.utils import logger


def convert_raw_files(
        raw_file_paths: List[Path],
        output_folder: Optional[Path] = None,
        num_threads=1,
        ms_level: str = "2-"):
    """
    Converts raw files to mzML files using ThermoRawFileParser.
    :param raw_file_paths: List of paths to raw files
    :param output_folder: Path to output folder
    :param num_threads: Number of threads to use
    :param ms_level: MS level to convert; keep at default if no specific MS level is needed for further processing
    :return: None
    """
    convert_raw_mzml_batch(raw_files=raw_file_paths, output_folder=output_folder, num_threads=num_threads,
                           ms_level=ms_level)


def convert_and_get_path(raw_type, threads, raw_dir, msms, output_dir):
    """
    Convert raw files to mzML and return the path to the mzML directory
    :param raw_type: 'raw' or 'mzml'
    :param threads: Number of threads to use
    :param raw_dir: Directory containing input files
    :param msms: msms.txt file; used to get a list of all involved raw files
    :param output_dir: Directory to generate mzML folder in and save files after conversion
    :return: Path to mzML directory
    """
    if raw_type == 'raw':
        raw_file_paths = [raw_dir / f'{f}.raw' for f in list(set(msms['Raw file']))]
        logger.info(f'Converting .raw files to mzML')
        if output_dir / 'mzml' not in output_dir.iterdir():
            (output_dir / 'mzml').mkdir()
        convert_raw_files(raw_file_paths, output_dir / 'mzml', threads)
        mzml_dir = output_dir / 'mzml'
    elif raw_type == 'mzml':
        mzml_dir = raw_dir
    else:
        raise ValueError(f'Unknown raw type: {raw_type}')
    return mzml_dir
