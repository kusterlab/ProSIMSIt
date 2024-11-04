import pandas as pd
from pathlib import Path
from typing import List

from oktoberfest import preprocessing as pp


SIMSI_COLUMNS = []


def read_msms_singlecol(msms_path: Path, onlycolumn):
    if msms_path.is_dir():
        msms_path = msms_path / 'msms.txt'
    return pd.read_csv(msms_path, sep='\t', usecols=[onlycolumn])


def read_fasta(fasta_path):
    with open(fasta_path) as file:
        return file.read()


def export_intermediate_file(df: pd.DataFrame, filename: str, location: Path):
    location = location / filename
    df.to_csv(location, sep='\t', index=False)


def generate_spectra_list(mzml_list: List[Path]):
    spectra_list = []
    for spectra_file in mzml_list:
        spectra = pp.load_spectra(filenames=spectra_file, parser="pyteomics")
        spectra_list.append(spectra)
    return spectra_list