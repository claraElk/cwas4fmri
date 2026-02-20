import os
import sys
import argparse
import pandas as pd
import numpy as np

from pathlib import Path
from .worklfow import workflow

def global_parser():
    """Create and return argument parser for CWAS analysis."""
    parser = argparse.ArgumentParser(description="Run CWAS Analysis")

    # BIDS app required arguments
    parser.add_argument(
        "bids_dir",
        action="store",
        type=Path,
        help="The directory with the input dataset (e.g. fMRIPrep derivative)formatted according to the BIDS standard",
    )
    parser.add_argument(
        "output_dir",
        action="store",
        type=Path,
        help="The directory where the output files should be stored",
    )
    parser.add_argument(
        "analysis_level",
        help="Level of the analysis that will be performed. Only group level is available",
        choices=["group"],
    )
    parser.add_argument(
        "--verbosity",
        help="""
        Verbosity level.
        """,
        required=False,
        choices=[0, 1, 2, 3],
        default=2,
        type=int,
        nargs=1,
    )
    parser.add_argument("--strategy", type=str, required=True, help="Strategy name")
    parser.add_argument("--phenotype", type=str, required=True, help="Path to phenotype file after QC")
    parser.add_argument("--atlas", type=str, required=True, help="Atlas name")
    parser.add_argument("--atlas_file", type=str, required=True, help="Path to ROI labels file")
    parser.add_argument("--patient", type=str, required=True, help="Name of the case group")
    parser.add_argument("--control", type=str, required=True, help="Name of the control group")
    parser.add_argument("--categorical_covariates", type=str, required=False, help="Columns with categorical covariates to include in the model, separated by commas")
    parser.add_argument("--numerical_covariates", type=str, required=False, help="Columns with numerical covariates to include in the model, separated by commas")
    
    return parser


def main(argv=None):
    args = global_parser().parse_args(argv)
    workflow(args)

if __name__ == "__main__":
    main(sys.argv[1:])