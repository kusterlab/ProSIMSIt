PROSIT_CONFIG = {
    "type": "Rescoring",
    "inputs": {
        "search_results_type": "Maxquant",
        "spectra_type": "mzml",
    },
    "fdr_estimation_method": "percolator",
    "allFeatures": False,
    "regressionMethod": "spline",
    "ssl": False,
    "ce_alignment_options": {
        "ce_range": [19, 50],
        "use_ransac_model": False
    }
}
