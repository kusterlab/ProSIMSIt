import logging

import pandas as pd

from picked_group_fdr import picked_group_fdr
from picked_group_fdr.pipeline import update_evidence_from_pout

logger = logging.getLogger(__package__ + "." + __file__)

def run_picked_protein_group_fdr(percolator_dir, picked_dir, fasta, enzyme):
    def add_extra_dashes_for_percolator(path_to_file, path_to_file_out):
        """
        deletes extra dashes from the epercolator output
        :param path_to_file:
        :param path_to_file_out:
        :return:
        """
        df = pd.read_csv(path_to_file, sep='\t')
        df['PSMId'] += '-1'
        df.to_csv(path_to_file_out, sep='\t', index=False)

    add_extra_dashes_for_percolator(f'{percolator_dir}/rescore_all.percolator.psms.txt',
                                    f'{percolator_dir}/rescore_all.percolator.psms.dash.txt')
    add_extra_dashes_for_percolator(f'{percolator_dir}/rescore_all.percolator.decoy.psms.txt',
                                    f'{percolator_dir}/rescore_all.percolator.decoy.psms.dash.txt')

    update_evidence_from_pout.main([
        '--mq_evidence', f'{picked_dir}/evidence.txt',
        '--perc_results', f'{percolator_dir}/rescore_all.percolator.psms.dash.txt', f'{percolator_dir}/rescore_all.percolator.decoy.psms.dash.txt',
        '--mq_evidence_out', f'{picked_dir}/updated_evidence.txt',
        '--pout_input_type', 'prosit'])


    if type(fasta) == list:
        picked_group_fdr.main([
            '--mq_evidence', f'{picked_dir}/updated_evidence.txt',
            '--protein_groups_out', f'{picked_dir}/group_results.txt',
            '--fasta', *fasta,
            '--methods', 'picked_protein_group_mq_input',
            '--enzyme', enzyme])
    else:
        picked_group_fdr.main([
            '--mq_evidence', f'{picked_dir}/updated_evidence.txt',
            '--protein_groups_out', f'{picked_dir}/group_results.txt',
            '--fasta', fasta,
            '--methods', 'picked_protein_group_mq_input',
            '--enzyme', enzyme])