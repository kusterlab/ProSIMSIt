import pandas as pd
from pathlib import Path


def read_msms_singlecol(msms_path: Path, onlycolumn):
    if msms_path.is_dir():
        msms_path = msms_path / 'msms.txt'
    return pd.read_csv(msms_path, sep='\t', usecols=[onlycolumn])