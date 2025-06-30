import numpy as np
import pandas as pd

import glob
import os
import json
from collections import defaultdict
from tqdm import tqdm

from .files import report_file

def filter_by_qc(json_file_path, pheno_filtered, out_p):
    """
    Process ratings and clean phenotype table by removing subjects with bad ratings
    """
    print("⏳ Identify subjects with bad QC ...")

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
    json_path = os.path.join(out_p, 'cwas_report.json')
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

    print("\n=== Summary QC ===")
    for key, value in summary_data.items():
        if isinstance(value, list):
            print(f"{key}: {len(value)}")
        else:
            print(f"{key}: {value}")
            
    print(f"✅ Information saved in:", json_path)

    return pheno_filtered_qc

def filter_by_fd(pheno_filtered_qc, derivatives_p, confounds_json, out_p, session, task, run, feature):
    """
    Filter subjects based on framewise displacement (FD)
    """
    print("\n⏳ Reject subjects based on mean FD>0.5 ...")
    print("This might take a moment, please do not interupt the process ...\n")
    
    # Find subjects processed by HALFpipe and collect their FD values
    mean_fd_list = []
    
    for _, row in tqdm(pheno_filtered_qc.iterrows()):
        json_file_path = os.path.join(derivatives_p, row['participant_id'], 'ses-{}'.format(session), "func",
                                                confounds_json.format(row['participant_id'], session, task, run, feature))
        if os.path.exists(json_file_path) :
            # Read the JSON file
            with open(json_file_path, 'r') as file:
                data = json.load(file)  

            mean_fd_list.append(data["MeanFramewiseDisplacement"])
        
    # Add mean FD values to phenotype dataframe
    pheno_filtered_qc.loc[:, 'mean_fd'] = mean_fd_list

    # Filter out subjects with high mean framewise displacement (FD >= 0.5)
    pheno_filtered_fd_mean = pheno_filtered_qc[pheno_filtered_qc['mean_fd'] < 0.5]
    subjects_with_mean_rejection = pheno_filtered_qc[pheno_filtered_qc['mean_fd'] > 0.5]['participant_id']
    
    # Define summary data
    summary_data = {
        "Total subjects before FD rejection": len(pheno_filtered_qc),
        "N Subjects with mean FD>0.5": len(subjects_with_mean_rejection),
        "Subjects with mean FD>0.5": sorted(subjects_with_mean_rejection),
    }

    # Save summary to file
    json_path = os.path.join(out_p, 'cwas_report.json')
    report_file(out_p, summary_data)

    print("\n=== Summary FD rejection ===")
    for key, value in summary_data.items():
        if isinstance(value, list):
            print(f"{key}: {len(value)}") 
        else:
            print(f"{key}: {value}")
            
    print(f"\n✅ Information saved in:", json_path)
    
    return pheno_filtered_fd_mean