import os
import argparse
import pandas as pd
import numpy as np

from pathlib import Path
from utils.tools import run_cwas_analysis

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run CWAS Analysis")
    parser.add_argument("--feature_settings", type=str, required=True, help="Path to feature settings file")
    parser.add_argument("--pheno_filtered_path", type=str, required=True, help="Path to filtered phenotype QC file")
    parser.add_argument("--derivatives_path", type=str, required=True, help="Path to derivatives data")
    parser.add_argument("--atlas_file", type=str, required=True, help="Path to ROI labels file")
    parser.add_argument("--out_p", type=str, required=True, help="Output path for results")
    parser.add_argument("--case_name", type=str, required=True, help="Name of the case group")
    parser.add_argument("--control_name", type=str, required=True, help="Name of the control group")
    parser.add_argument("--sequence_col", type=str, required=False, help="Column name for sequence data")
    parser.add_argument("--machine_col", type=str, required=False, help="Column name for machine data")

    args = parser.parse_args()
    print(f"Results will be saved to {args.out_p}.")

    
    # Ensure output directory exists
    os.makedirs(args.out_p, exist_ok=True)
    
    # Placeholder for CWAS analysis logic
    print(f"Running CWAS analysis on {args.feature_settings}")
    print(f"Results will be saved to {args.out_p}.")
    
    labels = pd.read_csv(args.atlas_file, sep='\t', header=None)
    roi_labels = labels[1].to_list()
    
    pheno_filtered_qc_fd = pd.read_csv(args.pheno_filtered_path, sep='\t', dtype={"participant_id": str})
    
    results = run_cwas_analysis(
        feature=args.feature_settings, 
        pheno_filtered_qc_fd=pheno_filtered_qc_fd, 
        derivatives_path=Path(args.derivatives_path), 
        roi_labels=roi_labels, 
        out_p=args.out_p, 
        case_name=args.case_name, 
        control_name=args.control_name,
        sequence_col=args.sequence_col, 
        machine_col=args.machine_col, 
        atlas="yeo2011"
    )