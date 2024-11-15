import logging
import glob
from pathlib import Path
from typing import List, Optional

import pandas as pd
import numpy as np

from simsi_transfer.thermo_raw import convert_raw_mzml_batch

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + __file__)


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


def prosit_to_simsi(path_to_msms, path_to_percolator, path_out, raw_file_hyphens=0):
    """
    Prepare a file usable as input for SIMSI-Transfer from the results of MaxQuant and Oktoberfest
    :param path_to_msms: Path to msms.txt file from MaxQuant
    :param path_to_percolator: Path to Oktoberfest output/results/percolator folder
    :param path_out: Path to save the output file
    :param raw_file_hyphens: Number of hyphens in the raw file name; required to properly split PSMId information
    :return: None
    """
    if path_out.is_file():
        logger.info(f'File {path_out} already exists; skipping new file generation')
        return

    prosit_target = pd.read_csv(path_to_percolator / 'rescore.percolator.psms.txt', sep="\t")
    prosit_decoy = pd.read_csv(path_to_percolator / 'rescore.percolator.decoy.psms.txt', sep="\t")

    prosit_all = pd.concat([prosit_target, prosit_decoy])
    prosit_all["Scan number"] = prosit_all['PSMId'].str.split('-').str[raw_file_hyphens + 1].astype(
        int)
    prosit_all["Raw file"] = prosit_all['PSMId'].str.split('-').str[0:raw_file_hyphens + 1]
    prosit_all["Raw file"] = prosit_all["Raw file"].apply(lambda raw_file: '-'.join(raw_file))
    prosit_all = prosit_all.loc[prosit_all["q-value"] <= 0.01]
    prosit_all = prosit_all[["Raw file", "Scan number", "posterior_error_prob", "score"]]

    msms100perc = pd.read_csv(path_to_msms, sep="\t")

    merged_df = prosit_all.merge(msms100perc, how='left', on=['Raw file', 'Scan number'], validate='1:1')
    merged_df["PEP"] = merged_df["posterior_error_prob"]
    merged_df["Score"] = merged_df["score"]
    merged_df = merged_df.drop("posterior_error_prob", axis=1)
    merged_df = merged_df.drop("score", axis=1)

    merged_df.to_csv(path_out, sep='\t', index=False)
    logger.info(f'Done preparing; saved SIMSI-ready file to {path_out}')


def prepare_input_for_second_oktoberfest(simsi_output):
    """
    Prepare a file for the second Oktoberfest run
    :param simsi_output: Path to SIMSI-Transfer output directory
    :return: Path to the msms.txt file usable as Oktoberfest input
    """
    msms_df = pd.read_csv(simsi_output / 'summaries/p10/p10_msms.txt', sep='\t')

    msms_df = msms_df[msms_df['identification'] == 't']
    msms_df = msms_df.rename(columns={'scanID': 'Scan number'})
    msms_df['Charge'] = msms_df['Charge'].astype(int)
    msms_df['Scan event number'] = msms_df['Scan number']
    msms_df.loc[msms_df['Score'].isna(), ['Score']] = 1
    msms_df['Mass'] = msms_df['Mass'].fillna((msms_df['m/z'] - 1.0078 + 0.0005) * msms_df['Charge'])

    msms_for_oktoberfest = simsi_output / 'summaries/p10/msms.txt'
    msms_df.to_csv(msms_for_oktoberfest, sep='\t', index=False)
    return msms_for_oktoberfest


def merge_rescore_files(ok1_dir, ok2_dir: Path, output_dir: Path):
    """
    Merge rescore files from the first and second Oktoberfest runs
    :param ok1_dir: Output directory of the first Oktoberfest run
    :param ok2_dir: Output directory of the second Oktoberfest run
    :param output_dir: Directory to save the merged file
    :return: None
    """
    rescoretab = pd.DataFrame()
    if (output_dir / 'rescore_all.tab').is_file():
        logger.info('rescore_all.tab already exists; reusing')
        return
    else:
        for f in glob.iglob(str(ok2_dir / '*rescore.tab')):
            temp = pd.read_csv(f, sep='\t')
            rescoretab = pd.concat([rescoretab, temp])
    temp = pd.read_csv(ok1_dir / 'rescore.tab', sep='\t').drop(columns=['ExpMass'])
    rescoretab = pd.concat([rescoretab, temp])
    rescoretab.insert(loc=4, column="ExpMass", value=rescoretab.groupby(["filename", "ScanNr"]).ngroup())
    rescoretab.to_csv(output_dir / 'rescore_all.tab', sep='\t', index=False)


def translate_modified_sequences_in_psmid(inpseries):
    """
    Translate modified sequences from MaxQuant to SIMSI format to be used in PSMId columns
    :param inpseries: Series containing PSMId's from Oktoberfest
    :return: Series containing PSMId's with modifications translated into MaxQuant format
    """
    split = inpseries.str.split('-')
    return (split.str[0] + '-' + split.str[1] + '-_' + split.str[3]
            .str.replace('[UNIMOD:737]', '', regex=False)
            .str.replace('[UNIMOD:35]', '(Oxidation (M))', regex=False)
            .str.replace('[UNIMOD:21]', '(Phospho (STY))', regex=False)
            .str.replace('[UNIMOD:4]', '', regex=False) + '_')


def prepare_for_building_evidence(path_to_percolator_result, path_to_percolator_decoy, path_to_simsi_msms,
                                  path_to_mq_msms100perc, path_to_mq_summary, path_to_output, number_of_hyphen=0):
    """
    Prepare a file in the shape of a simsi summary file, that includes all target and decoy PSMs generated during the workflow
    :param path_to_percolator_result: Path to the rescore.psms from the second Oktoberfest run
    :param path_to_percolator_decoy: Path to the rescore.decoy.psms from the second Oktoberfest run
    :param path_to_simsi_msms: Path to the msms.txt file generated by SIMSI-Transfer
    :param path_to_mq_msms100perc: Path to the 100% FDR msms.txt file from MaxQuant
    :param path_to_mq_summary: Path to the summary.txt file from MaxQuant
    :param path_to_output: Path to save the merged file
    :param number_of_hyphen: Number of hyphens in the raw file name; required to properly split PSMId information
    :return: None
    """
    if Path(path_to_output).is_file():
        logger.info(f'Merged msms.txt file found at {path_to_output}. Skipping file generation.')
        return

    all_PEPs = pd.DataFrame()
    deduplicated_PSMs = set()
    percolator = pd.read_csv(path_to_percolator_result, usecols=['PSMId', 'filename', 'posterior_error_prob'], sep='\t')
    percolator["Scan number"] = percolator['PSMId'].str.split('-').str[number_of_hyphen + 1].astype(int)
    percolator["ID"] = percolator['filename'] + '-' + percolator["Scan number"].astype(str)
    percolator["PSMId"] = translate_modified_sequences_in_psmid(percolator["PSMId"])
    deduplicated_PSMs = deduplicated_PSMs.union(set(percolator['PSMId']))
    all_ids = set(percolator['ID'])
    all_PEPs = pd.concat([all_PEPs, percolator[['ID', 'posterior_error_prob']]])
    del percolator
    percolator_decoys = pd.read_csv(path_to_percolator_decoy, usecols=['PSMId', 'filename', 'posterior_error_prob'],
                                    sep='\t')
    percolator_decoys["Scan number"] = percolator_decoys['PSMId'].str.split('-').str[number_of_hyphen + 1].astype(int)
    percolator_decoys["ID"] = percolator_decoys['filename'] + '-' + percolator_decoys["Scan number"].astype(str)
    percolator_decoys["PSMId"] = translate_modified_sequences_in_psmid(percolator_decoys["PSMId"])
    deduplicated_PSMs = deduplicated_PSMs.union(set(percolator_decoys['PSMId']))
    all_ids_decoys = set(percolator_decoys['ID'])
    all_PEPs = pd.concat([all_PEPs, percolator_decoys[['ID', 'posterior_error_prob']]])
    del percolator_decoys

    msms_simsi = pd.read_csv(path_to_simsi_msms, sep='\t')
    msms_simsi['ID'] = msms_simsi['Raw file'] + '-' + msms_simsi['scanID'].astype(str)
    msms_simsi['PSMId'] = msms_simsi['Raw file'] + '-' + msms_simsi['scanID'].astype(str) + '-' + msms_simsi[
        'Modified sequence']
    msms_simsi = msms_simsi[msms_simsi['PSMId'].isin(deduplicated_PSMs)]
    msms_simsi = msms_simsi.drop(columns=['PSMId'])
    keep_cols = set(msms_simsi.columns)
    keep_cols.add("Scan number")
    keep_cols = keep_cols - {'Fraction', 'MS scan number', 'clusterID', 'Experiment', 'mod_ambiguous', 'ID', 'PEP',
                             'summary_ID', 'identification', 'scanID', 'raw_ambiguous', 'Phospho (STY) Probabilities'}

    ids_not_in_simsi = all_ids - set(msms_simsi['ID'])
    decoys_not_in_simsi = all_ids_decoys - set(msms_simsi['ID'])

    msms100 = pd.read_csv(path_to_mq_msms100perc, usecols=list(keep_cols), sep="\t")
    msms100["ID"] = msms100["Raw file"] + '-' + msms100["Scan number"].astype(str)
    msms100 = msms100[msms100['ID'].isin(ids_not_in_simsi.union(decoys_not_in_simsi))]
    msms100 = msms100.rename(columns={"Scan number": "scanID"})

    summary = pd.read_csv(path_to_mq_summary, sep='\t')
    if 'Fraction' not in summary.columns:
        summary['Fraction'] = 1
    summary = summary[['Raw file', 'Experiment', 'Fraction']]

    msms100 = pd.merge(msms100, summary, on='Raw file', how='left')

    highest_summary = msms_simsi['summary_ID'].max()
    msms100['summary_ID'] = range(highest_summary + 1, highest_summary + 1 + len(
        msms100))  # has to be unique or build_evidence() doesn't work

    msms100['raw_ambiguous'] = np.nan
    msms100['mod_ambiguous'] = np.nan
    msms100['Phospho (STY) Probabilities'] = np.nan
    msms100['MS scan number'] = np.nan
    msms100['clusterID'] = np.nan
    msms100['identification'] = 'd'

    msms_simsi = pd.concat(
        [msms_simsi[msms_simsi['ID'].isin(all_ids) | msms_simsi['ID'].isin(all_ids_decoys)], msms100],
        ignore_index=True)
    msms_simsi = msms_simsi.merge(all_PEPs, on='ID', how='left', validate='1:1')
    msms_simsi.to_csv(path_to_output, sep='\t', index=False)


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
