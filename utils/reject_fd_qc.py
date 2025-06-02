import numpy as np
import pandas as pd

import glob
import os
import json
from collections import defaultdict
from tqdm import tqdm

def filter_by_qc(json_file_path, pheno_filtered, out_p):
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

def filter_by_fd(pheno_filtered_qc, derivatives_p, confounds_p, out_p):
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
    
    # Find subjects processed by HALFpipe and collect their FD values
    mean_fd_list = []
    max_fd_list = []
    
    for rid, row in tqdm(pheno_filtered_qc.iterrows()):
        confounds_path = glob.glob(os.path.join(derivatives_p, confounds_p.format(row['participant_id'], row['participant_id'])),
                                  recursive=True)
        
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