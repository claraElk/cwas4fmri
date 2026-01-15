# Sebastian Urchs 

#from statsmodels.sandbox.stats.multicomp import multipletests as stm
import sys
import glob
import os
import json
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy import io as sio
from statsmodels.stats.multitest import multipletests as stm

from .stats import *

def summarize_glm(glm_table, mask_2d, roi_labels):
    """
    Summarize GLM results:
    - Converts flattened upper-triangle edges back to full square matrix
    - Computes FDR q-values
    """

    out_table = glm_table.copy()

    # Compute qval (pvals FDR corrected)
    # Because of NaN we need to create a mask
    pvals = glm_table.pvals.values
    mask = np.isfinite(pvals)

    qval = np.full_like(pvals, np.nan, dtype=float)
    _, qval_valid, _, _ = stm(pvals[mask], alpha=0.05, method="fdr_bh")
    qval[mask] = qval_valid
    out_table["qval"] = qval

    # Convert flattened betas/pvals to full matrices
    beta_table = pd.DataFrame(conn2mat(out_table.betas.values, mask_2d), 
                                    index=roi_labels, columns=roi_labels)
    
    stand_beta_table = pd.DataFrame(conn2mat(out_table.stand_betas.values, mask_2d), 
                                    index=roi_labels, columns=roi_labels)
    
    pval_table = pd.DataFrame(conn2mat(out_table.pvals.values, mask_2d),
                              index=roi_labels, columns=roi_labels)
    
    qval_table = pd.DataFrame(conn2mat(out_table.qval.values, mask_2d),
                              index=roi_labels, columns=roi_labels)

    return out_table, stand_beta_table, qval_table, pval_table, beta_table

def average_runs(corr_files):
    """
    Average connectomes through runs
    Args:
        corr_files (list): list of path to 1 subject HALFpipe correlation matrices

    Returns:
        avg_mat (np): Averaged connectivity matrix of 1 subject (n roi x n roi)
    """
    
    # Load all matrices for this subject
    matrices = [pd.read_csv(cf, sep="\t", header=None, dtype=np.float32).values 
                for cf in corr_files]
    
    # Average across runs using nanmean
    avg_mat = np.nanmean(matrices, axis=0)

    return avg_mat

def fd_mean_extraction(json_files):
    """
    Extract the mean FD in HALFpipe json file and average through runs

    Args:
        json_files (list): list of path to 1 subject json file(s)

    Returns:
        np.float() : averaged mean_fd across runs
    """

    fd_values = []
    
    for jf in json_files:
        try:
            with open(jf, "r") as f:
                metadata = json.load(f)
            if "FDMean" in metadata:
                fd_values.append(metadata["FDMean"])
        except Exception as e:
            print(f"Warning: could not read JSON for {subj}: {e}")

    return np.mean(fd_values) if fd_values else np.nan

def process_connectivity_matrix(pheno_filtered_fd, feature, derivatives_path):
    """
    Process connectivity matrices and return CWAS-ready flattened matrices.
    - Averages multiple runs per subject
    - Keeps NaNs
    - Uses upper triangle without diagonal

    Args:
        pheno_filtered_fd (pd.DataFrame) : dataframe with only good QC subjects
            - Columns age, gender and diagnosis are mandatory
        feature (str): pipeline name
        derivative_path (str): path to HALFpipe derivatives
    """

    print("Processing connectivity matrices ...")

    flattened_matrices = []
    valid_subject_indices = []
    fdmean_values = []

    for idx, row in tqdm(pheno_filtered_fd.iterrows(), total=len(pheno_filtered_fd)):
        subj = str(row["participant_id"])
        
        if not subj.startswith("sub-"):
            # Check if participant_id starts with sub-
            subj = f"sub-{subj}"

        # Find all correlation matrices for this subject
        corr_pat = f"{subj}_*feature-{feature}_*desc-correlation_matrix.tsv"
        corr_files = list(Path(derivatives_path).rglob(corr_pat))
    
        # --- Average runs
        avg_mat = average_runs(corr_files)
        
        # Create upper triangle mask without diagonal
        n_rois = avg_mat.shape[0]
        mask_2d = np.triu(np.ones((n_rois, n_rois), dtype=bool), k=1)
        
        # Flatten upper triangle for CWAS
        flattened_matrices.append(avg_mat[mask_2d])
        valid_subject_indices.append(idx)

        # --- FDMean extraction ---
        json_pat = f"{subj}_*feature-{feature}_*_timeseries.json"
        
        json_files = list(Path(derivatives_path).rglob(json_pat))
        fdmean_values.append(fd_mean_extraction(json_files))

    # --- Build final conn_stack ---
    conn_stack = np.vstack(flattened_matrices)  # shape: (n_subjects, n_edges)
    pheno_filtered_fd = pheno_filtered_fd.loc[valid_subject_indices].reset_index(drop=True)
    # WARNING: this might fail, to be changed, should be using valid_subject_indices
    pheno_filtered_fd["mean_fd"] = fdmean_values

    print(f"Processed {len(pheno_filtered_fd)} subjects")
    print(f"Connectivity stack shape (n_subjects x n_edges): {conn_stack.shape}")

    return conn_stack, pheno_filtered_fd, mask_2d

def run_cwas_analysis(feature, 
                      pheno_filtered_qc_fd, 
                      derivatives_path, 
                      roi_labels, 
                      out_p, 
                      case_name, 
                      control_name,
                      sequence_col, 
                      machine_col, 
                      atlas="schaeferCombined", 
                     ):
    """
    Run CWAS analysis for all features defined in the JSON specification
    
    Args:
        feature (str): Pipeline name (as indicated in HALFpipe filename)
        pheno_filtered_qc_fd (pd.DataFrame): Filtered phenotype data (only good subjects)
        derivatives_p (str): Path to derivatives directory
        roi_labels (list): List of ROI labels
        out_p (str): Output directory path
        case_name (str): case variable as indicated in phenotype file
        control_name (str): control group variable as indicated in phenotype file
        sequence_col (str): if different sequence, indicate column name
        sequence_col (str): if different machine, indicate column name
        atlas (str): name of the atlas, as indicated in HALFpipe filename.
            - Default: schaeferCombined
    
    Returns:
        None
    """
    
    print("Run CWAS analysis for all features defined in the JSON specification ...")
    print("This will take a moment, please do not interupt the process ...")

    # Create output directory if it doesn't exist
    out_p = Path(out_p)
    out_p.mkdir(parents=True, exist_ok=True)
    
    # Define regressors
    # TODO: need to change that at some point to accept more variables
    # but works for now
    list_regressor = ['age', 'mean_fd', 'C(gender)']
    if sequence_col :
        list_regressor.append(f'C({sequence_col})')
    if machine_col  :
        list_regressor.append(f'C({machine_col})')
        
    regressors = ' + '.join(list_regressor)
    print(f'\nregressors used in the model: {regressors}')
    
    # Process each feature
    print(f"\nProcessing feature: {feature}")
    
    # Process connectivity matrix
    conn_stack, final_df, mask_2d = process_connectivity_matrix(
        pheno_filtered_fd=pheno_filtered_qc_fd,
        feature=feature,
        derivatives_path=derivatives_path,
    )
    
    # Perform CWAS analysis
    glm_con = glm_wrap_cc(conn_stack, 
                          final_df, 
                          group='diagnosis', 
                          case=case_name, 
                          control=control_name, 
                          regressors=regressors, 
                          report=True)
    
    # Get results
    table_con, table_stand_beta_con, table_qval_con, table_pval_con, table_beta_con = summarize_glm(
        glm_con, mask_2d, roi_labels
    )
    
    # Save results
    base_filename = f'cwas_{case_name}_{control_name}_rsfmri_{feature}_{atlas}'
    table_con.to_csv(out_p / f'{base_filename}.tsv', sep='\t')
    table_stand_beta_con.to_csv(out_p / f'{base_filename}_standardized_betas.tsv', sep='\t')
    table_qval_con.to_csv(out_p / f'{base_filename}_fdr_corrected_pvalues.tsv', sep='\t')
    table_pval_con.to_csv(out_p / f'{base_filename}_pvalues.tsv', sep='\t')
    table_beta_con.to_csv(out_p / f'{base_filename}_betas.tsv', sep='\t')
    
    print(f"Completed processing for feature: {feature}")
    print(f"\nAnalysis complete. Results saved to: {out_p}")