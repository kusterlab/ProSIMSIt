import os
import sys
import time
import subprocess
from datetime import datetime

from oktoberfest import runner
from simsi_transfer import main as simsi

from prosimsit.oktoberfest_functions import prepare_second_oktoberfest_run, preprocess_spectra_files, annotate_library, \
    generate_pred_files, calculate_featuers
from prosimsit.picked_fdr_functions import run_picked_protein_group_fdr
from . import __version__, __copyright__
from .utils import *
from .command_line_interface import read_config
from .simsi_functions import prepare_simsi_files, build_evidence
from .io import read_msms_singlecol

#logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def main(argv):
    config = read_config(argv)
    print(config)

    output_dir = Path(config['general']['output'])
    threads = config['general']['threads']

    maxquant_dir = Path(config['inputs']['maxquant_results'])
    raw_dir = Path(config['inputs']['spectra'])
    raw_type = config['inputs']['spectra_type']

    output_dir.mkdir(parents=True, exist_ok=True)

    module_name = ".".join(__name__.split(".")[:-1])
    file_logger = logging.FileHandler(output_dir / Path('ProSIMSIt.log'))
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    file_logger.setFormatter(formatter)
    logging.getLogger(module_name).addHandler(file_logger)

    starttime = datetime.now()

    logger.info(f'ProSIMSIt version {__version__}')
    logger.info(f'{__copyright__}')
    logger.info(f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}')

    logger.info(f'Starting ProSIMSIt')
    logger.info('')


    logger.info(f'Retrieving .raw files')
    msms = read_msms_singlecol(maxquant_dir, 'Raw file')

    mzml_dir = convert_to_mzml(raw_type, threads, raw_dir, msms, output_dir)

    logger.info(f'Building config.json for first Oktoberfest run')
    oktoberfest_config_path = output_dir / 'config_oktoberfest.json'
    generate_oktoberfest_config(config, mzml_dir, oktoberfest_config_path)

    logger.info(f'Executing first Oktoberfest run')
    if (output_dir / 'oktoberfest_1_out/results/percolator/rescore.percolator.psms.txt').is_file():
        logger.info(f'Found previous Oktoberfest run; skipping...')
    else:
        runner.run_job(oktoberfest_config_path)

    logger.info(f'Preparing input file for SIMSI-Transfer')
    raw_file_hyphen = os.listdir(mzml_dir)[0].count('-')
    os.makedirs(output_dir / 'simsi_input', exist_ok=True)

    ok1_percolator = output_dir / 'oktoberfest_1_out' / 'results' / 'percolator'

    prosit_to_simsi(
        maxquant_dir / 'msms.txt',
        ok1_percolator,
        output_dir / 'simsi_input' / 'msms.txt',
        raw_file_hyphens=raw_file_hyphen)

    prepare_simsi_files(output_dir, maxquant_dir)

    logger.info(f'Starting SIMSI-Transfer')

    simsi_stringency = config['simsi']['stringency']
    simsi_max_pep = config['simsi']['max_pep']
    tmt_level = config['general']['tmt_ms_level']
    simsi_output = output_dir / 'simsi_output'

    simsi_args = [
        '--mq_txt_folder', str(output_dir / 'simsi_input'),
        '--raw_folder', str(mzml_dir),
        '--output_folder', str(simsi_output),
        '--cache_folder', str(output_dir / 'simsi_output'),
        '--stringencies', str(simsi_stringency),
        '--maximum_pep', str(simsi_max_pep),
        '--num_threads', str(threads),
        '--tmt_ms_level', str(tmt_level),
        '--ambiguity_decision', 'keep_all',
        '--skip_evidence', '--skip_msmsscans'
    ]
    if (simsi_output / 'summaries/p10/p10_msms.txt').is_file():
        logger.info(f'Found previous SIMSI-Transfer run; skipping...')
    else:
        simsi.main(simsi_args)
    logger.info(f'Finished SIMSI-Transfer!')

    logger.info(f'Starting second Oktoberfest run')
    msms_for_prosit_2 = prepare_for_second_prosit_run(simsi_output)

    simsi_output = output_dir / 'simsi_output'
    dir_msms = msms_for_prosit_2

    conf = prepare_second_oktoberfest_run(mzml_dir, oktoberfest_config_path, dir_msms, output_dir)
    spectra_files = preprocess_spectra_files(conf)
    annotate_library(spectra_files, conf)
    generate_pred_files(conf)
    calculate_featuers(spectra_files, conf)
    logger.info(f'Finished second Oktoberfest run')

    #####

    logger.info(f'Preparing for percolator run')
    percolator_dir = output_dir / 'ProSIMSIt/percolator'
    os.makedirs(percolator_dir, exist_ok=True)
    merge_rescore_files(ok1_dir=ok1_percolator, ok2_dir=conf.output / 'results/percolator', output_dir=percolator_dir)

    #####

    logger.info(f'Starting Percolator run')

    target_psms = Path(f"{percolator_dir}/rescore_all.percolator.psms.txt")
    decoy_psms = Path(f"{percolator_dir}/rescore_all.percolator.decoy.psms.txt")
    target_peptides = Path(f"{percolator_dir}/rescore_all.percolator.peptides.txt")
    decoy_peptides = Path(f"{percolator_dir}/rescore_all.percolator.decoy.peptides.txt")
    input_file = Path(f"{percolator_dir}/rescore_all.tab")
    log_file = Path(f"{percolator_dir}/rescore_all.log")

    if target_psms.is_file():
        logger.info(f'Percolator run already exists; reusing')
    else:
        cmd = f"percolator --init-weights {ok1_percolator}/rescore.percolator.weights.csv \
                            --static \
                            --num-threads {threads} \
                            --subset-max-train 500000 \
                            --post-processing-tdc \
                            --testFDR 0.01 \
                            --trainFDR 0.01 \
                            --results-psms {target_psms} \
                            --decoy-results-psms {decoy_psms} \
                            --results-peptides {target_peptides} \
                            --decoy-results-peptides {decoy_peptides} \
                            {input_file} 2> {log_file}"

        subprocess.run(cmd, shell=True, check=True)
    logger.info(f'Finished Percolator run')

    #####

    config['general']['plotting'] = False
    if config['general']['plotting']:
        logger.info(f'Plotting summary plots...')
        # plotting
        logger.info(f'Finished plotting summary plots.')

    logger.info(f'Assembling evidence file for Picked Protein Group FDR')
    picked_dir = output_dir / 'ProSIMSIt/PickedProteinGroupFDR'
    os.makedirs(picked_dir, exist_ok=True)
    prepare_for_building_evidence(
        f'{percolator_dir}/rescore_all.percolator.psms.txt',
        f'{percolator_dir}/rescore_all.percolator.decoy.psms.txt',
        f'{simsi_output}/summaries/p10/p10_msms.txt',
        f'{maxquant_dir}/msms.txt',
        f'{maxquant_dir}/summary.txt',
        f'{picked_dir}/merged_msms.txt',
        number_of_hyphen=raw_file_hyphen)

    build_evidence(f'{picked_dir}/merged_msms.txt', maxquant_dir, picked_dir)
    logger.info(f'Evidence assembly finished!')

    logger.info(f'Applying Picked Protein Group FDR')
    run_picked_protein_group_fdr(percolator_dir, picked_dir, config['picked_protein_group_fdr']['fasta'], config['picked_protein_group_fdr']['enzyme'])
    # picked
    logger.info(f'Picked Protein Group FDR application finished!')

    endtime = datetime.now()
    logger.info(f'ProSIMSIt finished in {endtime - starttime}!')

    # TODO:
    # Filechecks to remove unnecessary re-writes


if __name__ == '__main__':
    main(sys.argv[1:])
