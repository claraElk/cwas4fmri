import os
import numpy as np
import pandas as pd
from pathlib import Path

def find_halfpipe_output(working_directory):
    reports_dir = os.path.join(working_directory, "reports") # Path to report folder
    subject_p = os.path.join(working_directory, "subject-list.txt") # Subject list

    derivatives_p = os.path.join(working_directory, "derivatives") # Path to derivative folder
    connectome_p = os.path.join(derivatives_p, 'halfpipe') # Path to halfpipe derivative folder

    json_spec_path = os.path.join(working_directory, "spec.json")  # Path to spec.json file
    json_exclude_qc_path = os.path.join(reports_dir, 'exclude.json') # Path to exclude.json file

    connectome_t = os.path.join('{}', '**', '*{}_*_feature-{}_atlas-{}_desc-correlation_matrix.tsv')
    confounds_p = os.path.join('fmriprep', '{}', '**', '{}_*_desc-confounds_timeseries.tsv')

    dict_halfpipe = {
        "reports_dir": reports_dir,
        "subject_p": subject_p,
        "json_spec_path": json_spec_path,
        "json_exclude_qc_path": json_exclude_qc_path,
        "derivatives_p": derivatives_p,
        "connectome_p": connectome_p,
        "connectome_t": connectome_t,
        "confounds_p": confounds_p}
    
    return dict_halfpipe
    
def verify_atals_files(atlas_file) :
    # 1. Verify atlas path
    print("Verifying altas location ...")
    print("path to access atlas:", atlas_file)
    
    if not os.path.exists(atlas_file):
        raise FileNotFoundError(f"Missing required file: {atlas_file}")
    
    try:
        # Verify file can be read as TSV
        with open(atlas_file, 'r') as f:
            first_line = f.readline()
            if not first_line:
                raise ValueError(f"Empty file: {atlas_file}")
    except Exception as e:
        raise ValueError(f"Error reading atlas file: {str(e)}")

    labels = pd.read_csv(atlas_file, sep='\t', header=None)
    conn_mask = np.tril(np.ones((len(labels),len(labels)))).astype(bool)
    roi_labels = labels[1].to_list()

    return labels, conn_mask, roi_labels
    

def verify_working_directory(working_directory) :
    # 2. Verify working directory
    print("\nVerifying working directory location ...")
    print("Path to access working directory:", working_directory)
    
    if not os.path.exists(working_directory):
        raise FileNotFoundError(f"Missing exclude.json in reports directory: {working_directory}")
    

def verify_exclude_json(json_exclude_qc_path, reports_dir) :
    # 3. Verify exclude.json in /reports
    print("\nVerifying exlude.json file location ...")
    print("Path to access reports folder:", reports_dir)

    if not os.path.exists(reports_dir):
        raise FileNotFoundError(f"Missing reports directory: {reports_dir}")
    
    exclude_file = json_exclude_qc_path
    print("Path to access exclude.json file:", exclude_file)
    
    if not os.path.exists(exclude_file):
        raise FileNotFoundError(f"Missing exclude.json in reports directory: {exclude_file}")
    
    print("\nAll required files and folder verified successfully")

    return True

def verify_files_and_directories(working_directory, json_exclude_qc_path, reports_dir):
    """
    Verify all required files and directories exist.
    
    Args:
        atlas_file (str): Path to the atlas file.
        working_directory (str): Path to the working directory.
        json_exclude_qc_path (str): Path to the exclude.json file.
    
    Returns:
        bool: True if all verifications pass, raises an error otherwise.
    """
    verify_working_directory(working_directory)
    verify_exclude_json(json_exclude_qc_path, reports_dir)
    
    return
