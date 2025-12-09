# Sebastian Urchs 

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

import json
import pandas as pd
import numpy as np
from tqdm import tqdm

# From other files
# from params import *
from .stats import *


def summarize_glm(glm_table, conn_mask, roi_labels):
    out_table = glm_table.copy()
    (fdr_pass, qval, _, _) = stm(glm_table.pvals, alpha=0.05, method='fdr_bh')
    out_table['qval'] = qval
    # Return to matrix form
    stand_beta_table = pd.DataFrame(conn2mat(out_table.stand_betas.values, conn_mask) , index=roi_labels, columns=roi_labels)
    qval_table = pd.DataFrame(conn2mat(out_table.pvals.values, conn_mask), index=roi_labels, columns=roi_labels)
    return out_table, stand_beta_table, qval_table

def extract_feature_settings(json_file_path):
    """
    Extract feature names and their settings from JSON file and create a table
    
    Args:
        json_file_path (str): Path to the JSON file containing feature settings
    
    Returns:
        tuple: (pd.DataFrame, str) containing feature settings and the common atlas name
    
    Raises:
        ValueError: If different atlases are found across features
    """
    print("Identify variables based on spec.json ...")
    # Read the JSON file
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    
    # Extract features and settings
    feature_settings = []
    feature_atlases = []
    
    for feature in data['features']:
        if feature['type'] == 'atlas_based_connectivity':
            feature_settings.append(feature['name'])
            feature_atlases.append(feature['atlases'][0])
            
    # Check if all features use the same atlas
    if len(set(feature_atlases)) > 1:
        raise ValueError(f"Error: Multiple different atlases found across features:\n{set(feature_atlases)}")
    
    # Get the common atlas name
    common_atlas = feature_atlases[0] if feature_atlases else None
    
    # Print results
    print("\nFeature Settings Found in spec.json:")
    print(feature_settings)
    
    print("\nAtlas Information:")
    if common_atlas:
        print(f"Common atlas used across all features: {common_atlas}")
    else:
        print("No atlas-based connectivity features found")
    
    return feature_settings, common_atlas

def process_ratings_and_clean_phenotype(json_file_path, pheno_filtered):
    """
    Process ratings and clean phenotype table by removing subjects with bad ratings
    
    Args:
        json_file_path (str): Path to the JSON file containing ratings
        phenotype_path (str): Path to the phenotype CSV file
        output_path (str): Path where to save the cleaned phenotype table
    """
    print("Identify subjects with bad QC ...")

    # Dictionary to store subject ratings
    subject_ratings = defaultdict(lambda: {'good': 0, 'bad': 0, 'uncertain': 0})
    
    # Read the JSON file
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    
    # Count ratings
    for entry in data:
        subject = entry.get('sub', 'Unknown')
        rating = entry.get('rating', 'unknown').lower()
        
        if rating == 'good':
            subject_ratings[subject]['good'] += 1
        elif rating == 'bad':
            subject_ratings[subject]['bad'] += 1
        elif rating == 'uncertain':
            subject_ratings[subject]['uncertain'] += 1
    
    # Create lists of subjects with different rating types
    subjects_with_bad = set()
    subjects_with_uncertain = set()
    subjects_without_bad = set()
    
    for subject, counts in subject_ratings.items():
        if counts['bad'] > 0:
            subjects_with_bad.add('sub-{}'.format(subject))
        elif counts['uncertain'] > 0:
            subjects_with_uncertain.add(subject)
        else:
            subjects_without_bad.add(subject)
    
    # Load phenotype table
    df = pheno_filtered.copy()
    
    # Filter out subjects with bad ratings
    pheno_filtered_qc = df[~df['participant_id'].isin(subjects_with_bad)]
    
    # Define summary data
    summary_data = {
        "Total subjects processed by HALFpipe": len(pheno_filtered),
        "N Subjects with bad ratinga": len(subjects_with_bad),
        "N Subjects with uncertain ratings": len(subjects_with_uncertain),
        "N Subjects with only good ratings": len(subjects_without_bad),
        "Total Subjects remaining after cleaning": len(pheno_filtered_qc),
        "Subjects with bad ratings": sorted(subjects_with_bad),

    }

    # Save summary to file
    # Load existing JSON if it exists
    json_path = os.path.join(out_p, 'cwas_subject_report.json')
    if os.path.exists(json_path):
        with open(json_path, 'r') as f:
            existing_data = json.load(f)
    else:
        existing_data = {}

    # Update existing data with new summary
    existing_data.update(summary_data)

    # Save updated JSON
    with open(json_path, 'w') as f:
        json.dump(existing_data, f, indent=4)

    # Optionally still print summary
    print("\n=== Summary QC ===")
    for key, value in summary_data.items():
        if isinstance(value, list):
            print(f"{key}: {len(value)}")  # Or print all items if preferred
        else:
            print(f"{key}: {value}")
            
    print(f"Information saved in:", json_path)

    return pheno_filtered_qc

def filter_subjects_by_fd(pheno_filtered_qc, derivatives_p):
    """
    Filter subjects based on framewise displacement (FD)
    
    Args:
        pheno_filtered_qc (pd.DataFrame): Filtered phenotype data
        derivatives_p (str): Path to derivatives directory
    
    Returns:
        tuple: (pheno_filtered_fd)
    """
    print("\nReject subjects based on mean FD>0.5 ...")
    print("This might take a moment, please do not interupt the process ...")

    connectome_t, confounds_p = check_path()
    
    # Find subjects processed by HALFpipe and collect their FD values
    mean_fd_list = []
    max_fd_list = []
    
    for rid, row in tqdm(pheno_filtered_qc.iterrows()):
        confounds_path = glob.glob(os.path.join(derivatives_p, confounds_p.format(row['participant_id'], row['participant_id'])))
        
        if any(os.path.exists(path) for path in confounds_path):
            max_fd = np.max(pd.read_csv(confounds_path[0], sep='\t')['framewise_displacement'])
            max_fd_list.append(max_fd)
            
            mean_fd = np.nanmean(pd.read_csv(confounds_path[0], sep='\t')['framewise_displacement'])
            mean_fd_list.append(mean_fd)
    
    # Add mean FD values to phenotype dataframe
    pheno_filtered_qc.loc[:, 'mean_fd'] = mean_fd_list
    pheno_filtered_qc.loc[:, 'max_fd'] = max_fd_list

    # Filter out subjects with high mean framewise displacement (FD >= 0.5)
    pheno_filtered_fd_mean = pheno_filtered_qc[pheno_filtered_qc['mean_fd'] < 0.5]
    subjects_with_mean_rejection = pheno_filtered_qc[pheno_filtered_qc['mean_fd'] > 0.5]['participant_id']
    
    
    print("\nReject subjects based on max FD>0.5 ...")
    
    # Filter out subjects with high max framewise displacement (FD >= 3.0)
    pheno_filtered_fd_mean_max = pheno_filtered_fd_mean[pheno_filtered_fd_mean['max_fd'] < 3.0]
    subjects_with_max_rejection = pheno_filtered_fd_mean[pheno_filtered_fd_mean['max_fd'] > 3.0]['participant_id']

    # Define summary data
    summary_data = {
        "Total subjects before FD rejection": len(pheno_filtered_qc),
        "N Subjects with mean FD>0.5": len(subjects_with_mean_rejection),
        "N Subjects with max FD>3.0": len(subjects_with_max_rejection),
        "Total Subjects remaining after FD cleaning": len(pheno_filtered_fd_mean_max),
        "Subjects with mean FD>0.5": sorted(subjects_with_mean_rejection),
        "Subjects with max FD>3.0": sorted(subjects_with_max_rejection),

    }

    # Save summary to file
    # Load existing JSON if it exists
    json_path = os.path.join(out_p, 'cwas_subject_report.json')
    if os.path.exists(json_path):
        with open(json_path, 'r') as f:
            existing_data = json.load(f)
    else:
        existing_data = {}

    # Update existing data with new summary
    existing_data.update(summary_data)

    # Save updated JSON
    with open(json_path, 'w') as f:
        json.dump(existing_data, f, indent=4)

    # Optionally still print summary
    print("\n=== Summary FD rejection ===")
    for key, value in summary_data.items():
        if isinstance(value, list):
            print(f"{key}: {len(value)}")  # Or print all items if preferred
        else:
            print(f"{key}: {value}")
            
    print(f"Information saved in:", json_path)
    
    return pheno_filtered_fd_mean_max


def process_connectivity_matrix(pheno_filtered_fd, feature, derivatives_path, conn_mask):
    """
    Process connectivity matrices for valid subjects and extract FDMean from timeseries JSON files.
    
    Args:
        pheno_filtered_fd (pd.DataFrame): Filtered phenotype data with FD < 0.5
        feature (str): Feature name
        derivatives_path (Path): Path to derivatives/halfpipe/
        conn_mask (array): Mask for connectome data
    
    Returns:
        tuple: (conn_stack, updated_pheno)
    """
    print("Process connectivity matrices for valid subjects ...")

    valid_subject_paths = []
    valid_subject_indices = []
    fdmean_values = []     # <--- store FDMean values

    for idx, row in tqdm(pheno_filtered_fd.iterrows(), total=len(pheno_filtered_fd)):
        subj = str(row["participant_id"])
    
        # Normalize: ensure it starts with "sub-"
        if not subj.startswith("sub-"):
            subj = f"sub-{subj}"
    
        # ---- Find all correlation matrices for this subject ----
        corr_pat = f"{subj}_*feature-{feature}_*desc-correlation_matrix.tsv"
        corr_files = list(derivatives_path.rglob(corr_pat))
    
        if not corr_files:
            fdmean_values.append(np.nan)
            continue
    
        # ---- Load and average correlation matrices ----
        matrices = []
        for cf in corr_files:
            mat = pd.read_csv(cf, sep="\t", index_col=0)
            matrices.append(mat)
    
        # Compute the average matrix
        if len(matrices) > 1:
            print("Average connectomes for", subj)
            avg_mat = sum(matrices) / len(matrices)
        else:
            avg_mat = matrices[0]
    
        # Save it somewhere (if needed)
        # You can store the averaged matrix so the rest of your pipeline uses only ONE
        # Here we just append the path of the first file, but you may want to change that
        valid_subject_paths.append(str(corr_files[0]))
        valid_subject_indices.append(idx)
    
        # ---- JSON metadata for FDMean ----
        json_pat = f"{subj}_*feature-{feature}_*_timeseries.json"
        json_files = list(derivatives_path.rglob(json_pat))
    
        # If multiple JSON files exist, average FDMean too
        fd_values = []
        for jf in json_files:
            try:
                with open(jf, "r") as f:
                    metadata = json.load(f)
                if "FDMean" in metadata:
                    fd_values.append(metadata["FDMean"])
            except Exception as e:
                print(f"Warning: could not read JSON for {subj}: {e}")
    
        if fd_values:
            fdmean_values.append(np.mean(fd_values))
        else:
            fdmean_values.append(np.nan)
            
    # Add FDMean column to phenotype table
    pheno_filtered_fd = pheno_filtered_fd.copy()
    pheno_filtered_fd["mean_fd"] = fdmean_values

    # ---- Stack connectomes ----
    conn_stack = np.array([
        pd.read_csv(p, sep="\t", header=None).values[conn_mask]
        for p in valid_subject_paths
    ])

    print(f"\nProcessing statistics for feature {feature}")
    print(f"Total subjects: {len(pheno_filtered_fd)}")
    print(f"Subjects with valid connectomes: {len(valid_subject_paths)}")

    # ---- Filter phenotype table so we keep ONLY subjects with connectomes ----
    pheno_filtered_fd = pheno_filtered_fd.loc[valid_subject_indices].reset_index(drop=True)
    pheno_filtered_fd["mean_fd"] = [
        fdmean_values[i] for i in valid_subject_indices
    ]

    return conn_stack, pheno_filtered_fd


def run_cwas_analysis(feature, pheno_filtered_qc_fd, 
                      derivatives_path, conn_mask, 
                      roi_labels, out_p, case_name, control_name,
                      sequence_col, medication_col, atlas="schaeferCombined"):
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
    
    
    # Define regressors
    list_regressor = ['age', 'C(gender)', 'mean_fd']
    if sequence_col :
        list_regressor.append(f'C({sequence_col})')
    if medication_col  :
        list_regressor.append('C(medication)')
        
    regressors = ' + '.join(list_regressor)
    print(f'\nregressors used in the model: {regressors}')
    
    # Process each feature
    print(f"\nProcessing feature: {feature}")
    
    # Process connectivity matrix
    conn_stack, final_df = process_connectivity_matrix(
        pheno_filtered_fd=pheno_filtered_qc_fd,
        feature=feature,
        derivatives_path=derivatives_path,
        conn_mask=conn_mask
    )

    # Rename case and control with 1 and o
    final_df = final_df.copy()
    final_df['diagnosis'] = final_df['diagnosis'].replace({case_name: 1, control_name: 0})
         
    # Perform CWAS analysis
    glm_con = glm_wrap_cc(conn_stack, final_df, 
                            group='diagnosis', case=1, control=0, 
                            regressors=regressors, report=True)
    
    # Get results
    table_con, table_stand_beta_con, table_qval_con = summarize_glm(
        glm_con, conn_mask, roi_labels
    )
    
    # Save results
    base_filename = f'cwas_{case_name}_{control_name}_rsfmri_{feature}_{atlas}'
    table_con.to_csv(out_p / f'{base_filename}.tsv', sep='\t')
    table_stand_beta_con.to_csv(out_p / f'{base_filename}_standardized_betas.tsv', sep='\t')
    table_qval_con.to_csv(out_p / f'{base_filename}_fdr_corrected_pvalues.tsv', sep='\t')
    
    print(f"Completed processing for feature: {feature}")

    print("\nAnalysis complete. Results saved to:")
    print(out_p)