import logging
import pandas as pd
import numpy as np
import psite_annotation as pa
from pathlib import Path
import re
from io import StringIO
import gc
from typing import Dict, List, Union, Tuple, Any, Optional
from pyascore import MassCorrector, SpectraParser, IdentificationParser, PyAscore

logger = logging.getLogger(__package__ + "." + __file__)


def preprocess_for_pyascore(df: pd.DataFrame) -> pd.DataFrame:
    """
    Preprocess peptide data for pyAscore analysis.

    Converts UNIMOD modifications from Percolator output to mass values and renames columns
    to match pyAscore expected format.

    :param df: DataFrame containing peptide data with columns:
            - peptide: peptide sequence with UNIMOD modifications
            - scan: scan number
            - Charge: peptide charge state
            - score: percolator score

    :return: DataFrame with columns: Raw file, scan, sequence, charge, percolator score
    """
    df = df.copy()

    # Define UNIMOD to mass replacements
    unimod_replacements = {
        '_.[UNIMOD:737]-': '',
        '._': '',
        '[UNIMOD:21]': '[79.966331]',  # Phosphorylation
        '[UNIMOD:737]': '',  # TMT removal; gets added as static modification later
        '[UNIMOD:4]': '[57.021464]',  # Carbamidomethylation
        '[UNIMOD:35]': '[15.9949]'  # Oxidation
    }

    # Apply replacements
    for unimod, replacement in unimod_replacements.items():
        df['peptide'] = df['peptide'].str.replace(unimod, replacement, regex=False)

    # Rename columns to match pyAscore expected format
    df = df.rename(columns={
        'peptide': 'sequence',
        'scan': 'scan',
        'Charge': 'charge',
        'score': 'percolator score'
    })

    df = df.sort_values(by=['Raw file', 'scan'])
    return df[['Raw file', 'scan', 'sequence', 'charge', 'percolator score']]


def pyascore_scoring_with_stringio(
        df: pd.DataFrame,
        mzml_path: Path,
        raw_file: str,
        current_counter: int,
        max_counter: int,
        mass_corrector: MassCorrector
) -> pd.DataFrame:
    """
    Run PyAscore analysis on a single raw file using StringIO for imitating disk IO operations.

    :param df: DataFrame containing PSM data for a single raw file
    :param mzml_path: Path to directory containing mzML files
    :param raw_file: Name of the raw file (without extension)
    :param current_counter: Current file counter for progress tracking
    :param max_counter: Total number of files to process
    :param mass_corrector: Configured pyAscore MassCorrector object

    :return: DataFrame containing PyAscore results with columns:
    Raw file, scan, localized_peptide, pepscore, ascores, alt_sites
    """

    print(f'Processing {raw_file}, {current_counter + 1} of {max_counter}')

    # Use StringIO for memory-efficient CSV processing
    csv_buffer = StringIO()
    try:
        df.to_csv(csv_buffer, index=False, sep='\t')
        csv_buffer.seek(0)
        id_parser = IdentificationParser(
            csv_buffer,
            'percolatorTXT',
            mass_corrector,
            static_mods={"n": 229.162932, "K": 229.162932, "C": 57.021464})
        psm_objects = id_parser.to_list()
    finally:
        csv_buffer.close()

    # Parse spectra file
    spectra_file = str(mzml_path / f'{raw_file}.mzML')
    spectra_parser = SpectraParser(spectra_file, 'mzML')
    spectra_objects = spectra_parser.to_dict()

    # Configure pyAscore for phosphorylation localization
    mod_mass = 79.966331  # Phosphorylation mass
    ascore = PyAscore(bin_size=100., n_top=10,
                      mod_group="STY",
                      mod_mass=mod_mass,
                      mz_error=.05,
                      fragment_types="by")
    ascore.add_neutral_loss('s', 97.9339)
    ascore.add_neutral_loss('t', 97.9339)

    pyascore_results = []
    for psm in psm_objects:
        # Check for phosphorylation modifications
        mod_select = np.isclose(psm["mod_masses"], mod_mass)
        nmods = np.sum(mod_select)

        if nmods >= 1:  # Process only phosphorylated peptides
            # Grab spectrum
            spectrum = spectra_objects[psm["scan"]]

            # Gather other (non-phospho) modifications
            aux_mod_pos = psm["mod_positions"][~mod_select].astype(np.uint32)
            aux_mod_masses = psm["mod_masses"][~mod_select].astype(np.float32)

            # Run PyAscore scoring
            ascore.score(mz_arr=spectrum["mz_values"],
                         int_arr=spectrum["intensity_values"],
                         peptide=psm["peptide"],
                         n_of_mod=np.sum(mod_select),
                         max_fragment_charge=max(1, psm["charge_state"] - 1),
                         aux_mod_pos=aux_mod_pos,
                         aux_mod_mass=aux_mod_masses)
            alt_sites = [",".join([str(site) for site in site_list]) for site_list in ascore.alt_sites]

            # Store scores for later use
            pyascore_results.append({
                "Raw file": raw_file,
                "scan": psm["scan"],
                "localized_peptide": ascore.best_sequence,
                "pepscore": ascore.best_score,
                "ascores": ";".join([str(s) for s in ascore.ascores]),
                "alt_sites": ";".join(alt_sites), })

    # Cleanup memory and garbage collection
    del ascore, spectra_parser, spectra_objects, id_parser, psm_objects
    gc.collect()

    return pd.DataFrame(pyascore_results)


def maxquantify_sequence(ascore_seq_column: pd.Series) -> pd.Series:
    """
    Convert PyAscore sequence format to MaxQuant-style modification notation.

    :param ascore_seq_column: Series containing pyAscore localized sequences

    :return: Series with MaxQuant-style modification notation
    """

    replacements = {
        r'n[229]': '',  # Remove N-terminal TMT
        r'n[458]': '',  # Remove N-terminal modifications
        r'[229]': '',  # Remove TMT
        r'[57]': '',  # Remove carbamidomethylation
        r'[80]': '(Phospho (STY))',  # Phosphorylation
        r'[16]': '(Oxidation (M))'  # Oxidation
    }

    result = ascore_seq_column.copy()
    for pattern, replacement in replacements.items():
        result = result.str.replace(pattern, replacement, regex=False)

    return result


def take_unique_else_raise(elements: List[Any]) -> Any:
    """
    Return the unique element from a list, raise error if multiple unique elements.

    :param elements: List of elements

    :return: The unique element if only one unique value exists

    :raises TypeError: If multiple unique elements are found
    """
    elements = set(elements)
    if len(elements) == 1:
        return list(elements)[0]
    else:
        raise TypeError(f'Multiple unique elements found; {elements}')


def csv_list_unique(elements: List[Any]) -> str:
    """

    Create a semicolon-separated string of unique elements.

    :param elements: List of elements to process

    :return: Semicolon-separated string of sorted unique elements
    """
    elements = list(set(elements))
    elements.sort()
    return ';'.join([str(i) for i in elements])


def process_ascores(series: pd.Series, function: str) -> str:
    """
    Process Ascore values using specified aggregation function.

    :param series: Series containing semicolon-separated Ascore values
    :param function: Aggregation function name ('mean', 'min', or 'max')

    :return: Semicolon-separated string of processed scores or NaN if no valid data

    :raises ValueError: If unsupported function is specified
    """

    valid_functions = {'mean', 'min', 'max'}
    if function not in valid_functions:
        raise ValueError(f"Unsupported function: {function}. Must be one of {valid_functions}")

    # Split each entry by semicolon and convert to float arrays
    score_arrays = []
    for entry in series:
        if pd.isna(entry):
            continue
        scores = [float(x) for x in str(entry).split(';')]
        score_arrays.append(scores)

    if not score_arrays:
        return np.nan

    score_matrix = np.array(score_arrays)

    # Apply aggregation function
    if function == 'mean':
        new_scores = np.mean(score_matrix, axis=0)
    elif function == 'min':
        new_scores = np.min(score_matrix, axis=0)
    elif function == 'max':
        new_scores = np.max(score_matrix, axis=0)

    return ';'.join(f"{score:.1f}" for score in new_scores)


def perform_pyascore(
    psms: pd.DataFrame,
    output_dir: Path,
    mzml_dir: Path
) -> pd.DataFrame:
    """

    Perform PyAscore analysis on all raw files in the dataset. Results are saved to 'ascore_table.txt' in the output
    directory.

    :param psms: DataFrame containing PSM data
    :param output_dir: Output directory path
    :param mzml_dir: Directory containing mzML files

    :return: DataFrame containing all PyAscore results
    """
    # Define modification masses
    modifications = {
        "n": 229.162932,    # n-term TMT6plex
        "M": 15.9949,       # Methionine oxidation
        "S": 79.966331,     # Serine phosphorylation
        "T": 79.966331,     # Threonine phosphorylation
        "Y": 79.966331,     # Tyrosine phosphorylation
        "C": 57.021464,     # Cysteine carbamidomethylation
        "K": 229.162932}    # Lysine TMT6plex
    mass_corrector = MassCorrector(modifications, mz_tol=1.5)
    preprocessed_df = preprocess_for_pyascore(psms)

    results = []
    counter = 0
    maxcounter = len(preprocessed_df['Raw file'].unique())

    logger.info(f'Processing {maxcounter} raw files with pyAscore')

    for raw_file, group in preprocessed_df.groupby('Raw file'):
        result = pyascore_scoring_with_stringio(group, mzml_dir, raw_file, counter, maxcounter, mass_corrector)
        results.append(result)
        counter += 1

    final_results = pd.concat(results, ignore_index=True)
    final_results.to_csv(output_dir / 'ProSIMSIt/ascore_table.txt', sep='\t', index=False)
    return final_results


def generate_phospho_table(
    output_dir: Path,
    mzml_dir: Path,
    fasta_path: Path
) -> None:
    """
    Generate a comprehensive phosphopeptide table with quantification and localization data.

    This function combines PSM data, quantification data, and PyAscore results to create a final phosphopeptide table
    with site localization information.

    :param output_dir: Output directory containing analysis results
    :param mzml_dir: Directory containing mzML files for PyAscore analysis
    :param fasta_path: Path to the protein FASTA file for site annotation
    """

    # Load PSM data
    psm_file = output_dir / 'ProSIMSIt/percolator/rescore_all.percolator.psms.txt'
    decoy_file = output_dir / 'ProSIMSIt/percolator/rescore_all.percolator.decoy.psms.txt'
    psms = pd.read_csv(psm_file, sep='\t')
    psms = pd.concat([psms, pd.read_csv(decoy_file, sep='\t')])

    # Load quan data
    quan_file = output_dir / 'ProSIMSIt/PickedProteinGroupFDR/merged_msms.txt'

    # First read the header to get column names, then read only relevant columns; needed to handle different TMT plexes
    quan_summary = pd.read_csv(quan_file, sep='\t', nrows=0)
    usecols = (
            [col for col in quan_summary.columns if 'Reporter intensity corrected' in col] +
            ['Raw file', 'scanID', 'Experiment', 'Fraction', 'Modified sequence']
    )

    quan_summary = pd.read_csv(quan_file, sep='\t', usecols=usecols)
    quan_summary['Fraction'] = quan_summary['Fraction'].astype('int16')

    # Process PSM data
    psms = psms.rename(columns={'filename': 'Raw file', 'proteinIds': 'Proteins'})
    psms = psms.loc[psms['q-value'] < 0.01]

    # Extract scan and charge information from PSMId
    psms['scan'] = psms['PSMId'].str.split('-').str[-4].astype('int32')
    psms['Charge'] = psms['PSMId'].str.split('-').str[-1].astype('int8')
    psms['Phosphorylations'] = psms['peptide'].str.count(r'\[UNIMOD\:21\]')

    # Run or load PyAscore results and merge with PSMs
    if (output_dir / 'ProSIMSIt/ascore_table.txt').is_file():
        logger.info(f'Existing pyAscore results found; reusing {output_dir / "ProSIMSIt/ascore_table.txt"}')
        final_results = pd.read_csv(output_dir / 'ProSIMSIt/ascore_table.txt', sep='\t')
    else:
        final_results = perform_pyascore(psms, output_dir, mzml_dir)

    psms = pd.merge(psms, final_results, on=['Raw file', 'scan'], how='left', validate='one_to_one')

    # Standardize reporter intensity column names (pad single digits with zero) to make them string-sortable and merge
    quan_summary = quan_summary.rename( columns={
        colname: re.sub(r' (\d)$', r' 0\1', colname)
        for colname in quan_summary.columns
        if 'Reporter intensity corrected' in colname
    }).rename(columns={'scanID': 'scan'})

    psms = pd.merge(psms, quan_summary, on=['Raw file', 'scan'], how='left', validate='one_to_one')

    # Process sequences
    psms = psms.rename(columns={'Modified sequence': 'MaxQuant sequence'})
    psms['Modified sequence'] = maxquantify_sequence(psms['localized_peptide'])

    psms['scanID'] = psms['Raw file'] + '|' + psms['scan'].astype(str)

    # Column aggregation for peptide table
    reporter_columns = [col for col in psms.columns if 'Reporter intensity corrected' in col]

    column_aggregation = {
        **{col: pd.NamedAgg(column=col, aggfunc='sum') for col in reporter_columns},
        **{'Number of PSMs': pd.NamedAgg(column='scanID', aggfunc='count')},
        **{'max ' + col: pd.NamedAgg(column=col, aggfunc='max') for col in ['pepscore']},
        **{'min ' + col: pd.NamedAgg(column=col, aggfunc='min') for col in ['pepscore']},
        **{'mean ' + col: pd.NamedAgg(column=col, aggfunc='mean') for col in ['pepscore']},
        **{'max ' + col: pd.NamedAgg(column=col, aggfunc=lambda x: process_ascores(x, 'max')) for col in ['ascores']},
        **{'min ' + col: pd.NamedAgg(column=col, aggfunc=lambda x: process_ascores(x, 'min')) for col in ['ascores']},
        **{'mean ' + col: pd.NamedAgg(column=col, aggfunc=lambda x: process_ascores(x, 'mean')) for col in ['ascores']},
        **{col + 's': pd.NamedAgg(column=col, aggfunc=csv_list_unique) for col in ['scanID', 'Fraction', 'Charge']},
    }

    # Create peptide-level summary
    peptide_table = psms.groupby(['Experiment', 'Modified sequence', 'Phosphorylations', 'Proteins']).agg(
        **column_aggregation).reset_index().sort_values(by=['Experiment', 'Modified sequence', 'Phosphorylations'])


    # Add peptide and phosphosite position annotations
    peptide_table = pa.addPeptideAndPsitePositions(peptide_table, fasta_path, pspInput=False)

    # Final sort and save
    peptide_table = peptide_table.sort_values(by=['Experiment', 'Modified sequence', 'Phosphorylations'])

    output_file = output_dir / 'ProSIMSIt/phosphopeptides.txt'
    peptide_table.to_csv(output_file, sep='\t', index=False)
    logger.info(f'Phosphopeptide table saved to {output_file}')
