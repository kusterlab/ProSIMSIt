import logging
import pandas as pd
import numpy as np
import psite_annotation as pa
from pathlib import Path
import re
from io import StringIO
import gc
from pyascore import MassCorrector, SpectraParser, IdentificationParser, PyAscore


logger = logging.getLogger(__package__ + "." + __file__)


def preprocess_for_pyascore(df):
    # Final columns:    scan sequence charge
    df = df.copy()
    # Replace UNIMOD with mod masses
    df['peptide'] = df['peptide'].str.replace('_.[UNIMOD:737]-', '', regex=False)
    df['peptide'] = df['peptide'].str.replace('._', '', regex=False)
    df['peptide'] = df['peptide'].str.replace('[UNIMOD:21]', '[79.966331]', regex=False)
    df['peptide'] = df['peptide'].str.replace('[UNIMOD:737]', '', regex=False)
    df['peptide'] = df['peptide'].str.replace('[UNIMOD:4]', '[57.021464]', regex=False)
    df['peptide'] = df['peptide'].str.replace('[UNIMOD:35]', '[15.9949]', regex=False)
    # rename columns
    df = df.rename(columns={'peptide': 'sequence', 'scan': 'scan', 'Charge': 'charge', 'score': 'percolator score'})
    df = df.sort_values(by=['Raw file', 'scan'])
    return df[['Raw file', 'scan', 'sequence', 'charge', 'percolator score']]


def pyascore_scoring_with_stringio(df, mzml_path, raw_file, current_counter, max_counter, mass_corrector):
    print(f'Processing {raw_file}, {current_counter + 1} of {max_counter}')
    csv_buffer = StringIO()
    try:
        df.to_csv(csv_buffer, index=False, sep='\t')
        csv_buffer.seek(0)
        id_parser = IdentificationParser(csv_buffer, 'percolatorTXT', mass_corrector,
                                         static_mods={"n": 229.162932, "K": 229.162932, "C": 57.021464})
        psm_objects = id_parser.to_list()
    finally:
        csv_buffer.close()

    spectra_file = str(mzml_path / f'{raw_file}.mzML')
    spectra_parser = SpectraParser(spectra_file, 'mzML')
    spectra_objects = spectra_parser.to_dict()

    mod_mass = 79.966331
    ascore = PyAscore(bin_size=100., n_top=10,
                      mod_group="STY",
                      mod_mass=mod_mass,
                      mz_error=.05,
                      fragment_types="by")
    ascore.add_neutral_loss('s', 97.9339)
    ascore.add_neutral_loss('t', 97.9339)

    pyascore_results = []
    for psm in psm_objects:
        # 4.1) Check for modification of interest
        mod_select = np.isclose(psm["mod_masses"], mod_mass)
        nmods = np.sum(mod_select)

        if nmods >= 1:
            # 4.2) Grab spectrum
            spectrum = spectra_objects[psm["scan"]]

            # 4.3) Gather other modifications into aux mods
            aux_mod_pos = psm["mod_positions"][~mod_select].astype(np.uint32)
            aux_mod_masses = psm["mod_masses"][~mod_select].astype(np.float32)

            # 4.4) Run scoring algorithm
            ascore.score(mz_arr=spectrum["mz_values"],
                         int_arr=spectrum["intensity_values"],
                         peptide=psm["peptide"],
                         n_of_mod=np.sum(mod_select),
                         max_fragment_charge=max(1, psm["charge_state"] - 1),
                         aux_mod_pos=aux_mod_pos,
                         aux_mod_mass=aux_mod_masses)
            alt_sites = [",".join([str(site) for site in site_list]) for site_list in ascore.alt_sites]

            # 4.5) Place scores into an object to use later
            pyascore_results.append({
                "Raw file": raw_file,
                "scan": psm["scan"],
                "localized_peptide": ascore.best_sequence,
                "pepscore": ascore.best_score,
                "ascores": ";".join([str(s) for s in ascore.ascores]),
                "alt_sites": ";".join(alt_sites), })

    # Clean up major objects
    del ascore, spectra_parser, spectra_objects, id_parser, psm_objects
    # Force garbage collection
    gc.collect()

    # 5) Make a dataframe of scores and write to a file
    # df = pd.merge(left=df, right=pd.DataFrame(pyascore_results), on=['Raw file', 'scan'], how='left', validate='one_to_one')
    return pd.DataFrame(pyascore_results)


def maxquantify_sequence(ascore_seq_column):
    mq_sequences = ascore_seq_column.str.replace(r'n[229]', '', regex=False).str.replace(r'n[458]', '',
                                                                                         regex=False).str.replace(
        r'[80]', '(Phospho (STY))', regex=False).str.replace(r'[16]', '(Oxidation (M))', regex=False).str.replace(
        r'[229]', '', regex=False).str.replace(r'[57]', '', regex=False)
    return mq_sequences


def take_unique_else_raise(elements):
    elements = set(elements)
    if len(elements) == 1:
        return list(elements)[0]
    else:
        raise TypeError(f'Multiple unique elements found; {elements}')


def csv_list_unique(elements):
    elements = list(set(elements))
    elements.sort()
    return ';'.join([str(i) for i in elements])


def process_ascores(series, function):
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
    if function == 'mean':
        new_scores = np.mean(score_matrix, axis=0)
    elif function == 'min':
        new_scores = np.min(score_matrix, axis=0)
    elif function == 'max':
        new_scores = np.max(score_matrix, axis=0)
    else:
        raise ValueError(f"Unsupported function: {function}")

    return ';'.join(f"{score:.1f}" for score in new_scores)


def perform_pyascore(psms, output_dir, mzml_dir):
    modifications = {"n": 229.162932,  # N-term TMT
                     "M": 15.9949,  # Methionine oxidation
                     "S": 79.966331,  # Serine Phoshorylation
                     "T": 79.966331,  # Threonine Phosphorylation
                     "Y": 79.966331,  # Tyrosine Phosphorylation
                     "C": 57.021464,  # Cysteine Carbamidomethylation
                     "K": 229.162932}  # Lysine TMT6plex
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
    final_results.to_csv(output_dir / 'ascore_table.txt', sep='\t', index=False)
    return final_results


def generate_phospho_table(output_dir, mzml_dir, fasta_path):
    psms = pd.read_csv(output_dir / 'ProSIMSIt/percolator' / 'rescore_all.percolator.psms.txt',
                       sep='\t')
    psms = pd.concat([psms, pd.read_csv(
        output_dir / 'ProSIMSIt/percolator' / 'rescore_all.percolator.decoy.psms.txt', sep='\t')])

    quan_summary = pd.read_csv(output_dir / 'PickedProteinGroupFDR/merged_msms.txt',
                               sep='\t', nrows=0)
    usecols = ([col for col in quan_summary.columns if 'Reporter intensity corrected' in col] +
               ['Raw file', 'scanID', 'Experiment', 'Fraction', 'Modified sequence'])
    quan_summary = pd.read_csv(output_dir / 'PickedProteinGroupFDR/merged_msms.txt',
                               sep='\t', usecols=usecols)


    quan_summary['Fraction'] = quan_summary['Fraction'].astype('int16')

    psms = psms.rename(columns={'filename': 'Raw file', 'proteinIds': 'Proteins'})
    psms = psms.loc[psms['q-value'] < 0.01]
    psms['scan'] = psms['PSMId'].str.split('-').str[-4].astype('int32')
    psms['Charge'] = psms['PSMId'].str.split('-').str[-1].astype('int8')
    psms['Phosphorylations'] = psms['peptide'].str.count(r'\[UNIMOD\:21\]')

    if (output_dir / 'ascore_table.txt').is_file():
        logger.info(f'Existing pyAscore results found; reusing {output_dir / "ascore_table.txt"}')
        final_results = pd.read_csv(output_dir / 'ascore_table.txt', sep='\t')
    else:
        final_results = perform_pyascore(psms, output_dir, mzml_dir)

    psms = pd.merge(psms, final_results, on=['Raw file', 'scan'], how='left', validate='one_to_one')

    quan_summary = quan_summary.rename(
        columns={colname: re.sub(r' (\d)$', r' 0\1', colname) for colname in quan_summary.columns if
                 'Reporter intensity corrected' in colname}).rename(columns={'scanID': 'scan'})

    psms = pd.merge(psms, quan_summary, on=['Raw file', 'scan'], how='left', validate='one_to_one')
    psms = psms.rename(columns={'Modified sequence': 'MaxQuant sequence'})
    psms['Modified sequence'] = maxquantify_sequence(psms['localized_peptide'])

    psms['scanID'] = psms['Raw file'] + '|' + psms['scan'].astype(str)

    column_aggregation = {
        **{col: pd.NamedAgg(column=col, aggfunc='sum') for col in psms.columns if 'Reporter intensity corrected' in col},
        **{'Number of PSMs': pd.NamedAgg(column='scanID', aggfunc='count')},
        **{'max ' + col: pd.NamedAgg(column=col, aggfunc='max') for col in ['pepscore']},
        **{'min ' + col: pd.NamedAgg(column=col, aggfunc='min') for col in ['pepscore']},
        **{'mean ' + col: pd.NamedAgg(column=col, aggfunc='mean') for col in ['pepscore']},
        **{'max ' + col: pd.NamedAgg(column=col, aggfunc=lambda x: process_ascores(x, 'max')) for col in ['ascores']},
        **{'min ' + col: pd.NamedAgg(column=col, aggfunc=lambda x: process_ascores(x, 'min')) for col in ['ascores']},
        **{'mean ' + col: pd.NamedAgg(column=col, aggfunc=lambda x: process_ascores(x, 'mean')) for col in ['ascores']},
        **{col + 's': pd.NamedAgg(column=col, aggfunc=csv_list_unique) for col in ['scanID', 'Fraction', 'Charge']},
    }

    peptide_table = psms.groupby(['Experiment', 'Modified sequence', 'Phosphorylations', 'Proteins']).agg(
        **column_aggregation).reset_index().sort_values(by=['Experiment', 'Modified sequence', 'Phosphorylations'])

    peptide_table = pa.addPeptideAndPsitePositions(
        peptide_table, fasta_path, pspInput=False)

    peptide_table = peptide_table.sort_values(by=['Experiment', 'Modified sequence', 'Phosphorylations'])
    peptide_table.to_csv(output_dir / 'phosphopeptides.txt', sep='\t', index=False)
