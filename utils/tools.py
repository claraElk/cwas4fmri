# TO BE DELETED

import numpy as np
import pandas as pd
from scipy import io as sio
from statsmodels.sandbox.stats.multicomp import multipletests as stm

# In python by default
import sys
import glob
import os
import json
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm

# From other files
from params import *
from .stats import *


def run_cwas_analysis(json_file_path, pheno_filtered_qc_fd, 
                      derivatives_p, connectome_p, conn_mask, 
                      roi_labels, out_p, case_name, control_name,
                      sequence_col, medication_col):
    """
    Run CWAS analysis for all features defined in the JSON specification
    
    Args:
        json_file_path (str): Path to JSON file containing feature specifications
        pheno_filtered_qc (pd.DataFrame): Filtered phenotype data
        derivatives_p (str): Path to derivatives directory
        connectome_p (str): Path to connectome directory
        conn_mask (array): Mask for connectome data
        roi_labels (list): List of ROI labels
        out_p (str): Output directory path
    
    Returns:
        None
    """
    print("Run CWAS analysis for all features defined in the JSON specification ...")
    print("This will take a moment, please do not interupt the process ...")

    # Create output directory if it doesn't exist
    out_p = Path(out_p)
    out_p.mkdir(parents=True, exist_ok=True)
    
    # Extract feature settings and validate atlases
    try:
        feature_settings, common_atlas = extract_feature_settings(json_file_path)
        
        if common_atlas:
            print(f"\nUsing common atlas: {common_atlas}")
        else:
            print("\nNo atlas-based connectivity features found")
            
    except ValueError as e:
        print(f"Error processing features: {str(e)}")
        sys.exit(1)    
    
    # Define regressors
    list_regressor = ['age', 'C(sex)', 'mean_fd', 'C(scanner)']
    if sequence_col :
        list_regressor.append('C(sequence)')
    if medication_col  :
        list_regressor.append('C(medication)')
        
    regressors = ' + '.join(list_regressor)
    print(f'\nregressors used in the model: {regressors}')
    
    # Process each feature
    for feature in tqdm(feature_settings):
        print(f"\nProcessing feature: {feature}")
        
        # Process connectivity matrix
        conn_stack, final_df = process_connectivity_matrix(
            pheno_filtered_fd=pheno_filtered_qc_fd,
            feature=feature,
            atlas=common_atlas,
            connectome_p=connectome_p,
            conn_mask=conn_mask
        )
        
        # Perform CWAS analysis
        glm_con = glm_wrap_cc(conn_stack, final_df, 
                             group='diagnosis', case=1, control=0, 
                             regressors=regressors, report=True)
        
        # Get results
        table_con, table_stand_beta_con, table_qval_con = summarize_glm(
            glm_con, conn_mask, roi_labels
        )
        
        # Save results
        base_filename = f'cwas_{case_name}_{control_name}_rsfmri_{feature}_{common_atlas}'
        table_con.to_csv(out_p / f'{base_filename}.tsv', sep='\t')
        table_stand_beta_con.to_csv(out_p / f'{base_filename}_standardized_betas.tsv', sep='\t')
        table_qval_con.to_csv(out_p / f'{base_filename}_fdr_corrected_pvalues.tsv', sep='\t')
        
        print(f"Completed processing for feature: {feature}")
    
    print("\nAnalysis complete. Results saved to:")
    print(out_p)