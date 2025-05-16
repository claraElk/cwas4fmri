# In python by default
import sys
import glob
import os
import json
import warnings
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm

# From other files
from .tools import *
from .stats import *
from params import *

# Required installation
import numpy as np
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)

def verify_table(phenotype_file_path, 
                 diagnosis_col, subject_col, 
                 age_col, sex_col, scanner_col,
                 sequence_col, medication_col,
                 case_name, control_name,
                 female_name, male_name):
    """
    Verify the phenotype file format and contents.
    
    Args:
        phenotype_file_path (str): Path to the phenotype file
        
    Returns:
        pandas.DataFrame: Verified phenotype file
        
    Raises:
        ValueError: If file format or content validation fails
    """
    # 1. Verify file extension
    print("Verifying location of phenotype file ...")
    print(phenotype_file_path)
    
    allowed_extensions = ['.tsv', '.csv', '.xlsx']
    file_extension = Path(phenotype_file_path).suffix.lower()
    
    if file_extension not in allowed_extensions:
        raise ValueError(f"Unsupported file extension: {file_extension}. "
                         f"Please use one of: {', '.join(allowed_extensions)}")
    
    # 2. Open the file according to extension
    try:
        if file_extension == '.tsv':
            df = pd.read_csv(phenotype_file_path, sep='\t')
        elif file_extension == '.csv':
            df = pd.read_csv(phenotype_file_path)
        else:  # .xlsx
            df = pd.read_excel(phenotype_file_path)
    except Exception as e:
        raise ValueError(f"Error reading file: {str(e)}")
    
    # 3. Verify column names
    print("\nVerifyng phenotype file content (column names, values) ...")

    col_label = {subject_col: "participant_id", 
                 diagnosis_col: "diagnosis", 
                 age_col: "age", 
                 sex_col: "sex",
                 scanner_col: 'scanner'
             }
    
    required_columns = [diagnosis_col, subject_col, age_col, sex_col]

    if sequence_col :
        col_label.update({sequence_col: 'sequence'})
        required_columns.append(sequence_col)

    if medication_col : 
        col_label.update({medication_col: 'medication'})
        required_columns.append(medication_col)

    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}. "
                         f"Required columns are: {', '.join(required_columns)}")
    
    # Rename columns
    print("Renaming columns ...")
    print(col_label)
    df.rename(columns=col_label, inplace=True)
    
    # Find scanner, sequence and medication
    scanner_name = np.unique(df['scanner'])
    print('Automatically found scanner IDs:', scanner_name)
    
    if sequence_col :
        sequence_name = np.unique(df['sequence'])
        print('Automatically found sequence IDs:', sequence_name)

    if medication_col : 
        medication_name = np.unique(df['medication'])
        print('Automatically found medications:', medication_name)

    
    # 4. Verify variables
    # Subject ID format
    with open(subject_p, "r") as f:
        subject_list = [line.strip() for line in f if line.strip()]  # Read and clean the lines

    # Create a regex pattern that matches any of the subjects
    subject_pattern = r'^sub-(?:' + '|'.join(subject_list) + r')$'

    invalid_ids = df[~df["participant_id"].str.match(subject_pattern)]

    if len(invalid_ids) > 0:
        warnings.warn(f"Invalid Subject ID format found. IDs should be in format 'sub-01' or 'sub-task01'. "
                      "Please verify that the following subjects were not included in HALFpipe subject-list.txt"
                      f"Invalid IDs: {invalid_ids['participant_id'].tolist()[:5]}")
    
    # Filter out rows with diagnosis values not matching the allowed categories
    allowed_diagnoses = [case_name, control_name]
    invalid_an = df[~df["diagnosis"].isin(allowed_diagnoses)]

    # Print and raise an error if any unexpected diagnosis values are found
    if not invalid_an.empty:
        print("WARNING: Invalid diagnosis values found:")
        print(invalid_an["diagnosis"].unique())  # Show only the unexpected diagnosis values
        
        # Change that as a WARNING
        warnings.warn(
            f"Phenotype file contains diagnoses outside of the expected set: {allowed_diagnoses}"
        )

    # Optionally, filter the dataframe to keep only valid rows
    df = df[df["diagnosis"].isin(allowed_diagnoses)]

    # Rename diagnosis value to be case:1 control:0
    label_case_control = {case_name : 0,
                         control_name : 1}
    
    print("Renaming diagnosis variables ...")
    print(label_case_control)
    
    df['diagnosis'] = df['diagnosis'].replace(label_case_control).infer_objects(copy=False)     
    # TODO: print count of case and controls
    
    # age verification (warning only)
    invalid_age = df[~df["age"].apply(lambda x: isinstance(x, (int, float)) and x > 0)]
    if len(invalid_age) > 0:
        print(f"Warning: Some age values may not be in years. "
              f"Invalid values found in rows: {invalid_age.index.tolist()[:5]}")
    
    # sex values
    invalid_sex = df[~df["sex"].isin([female_name, male_name])]
    if len(invalid_sex) > 0:
        raise ValueError(f"Invalid sex values found. Did not find {female_name} or {male_name} "
                         f"Please verify your sex variables")
    
    label_sex = {female_name : 0,
                 male_name : 1}
    
    print("Renaming sex variables ...")
    print(label_sex)
    
    df['sex'] = df['sex'].replace(label_sex).infer_objects(copy=False)     

    return df


def verify_excepted_files():
    """
    Verify the presence and format of required files in the working directory.
    
    Args:
        working_directory_folder (str): Path to the working directory
        
    Raises:
        FileNotFoundError: If required files are missing
        ValueError: If file contents are invalid
    """
    # 1. Verify atlas-Schaefer2018Combined_dseg.tsv
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
    
    # 2. Verify working directory
    print("\nVerifying working directory location ...")
    print("Path to access working directory:", working_directory)
    
    if not os.path.exists(working_directory):
        raise FileNotFoundError(f"Missing exclude.json in reports directory: {working_directory}")
    
    
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

def find_valid_subjects(pheno):
    confounds_t, confounds_p = check_path()
    
    print("Identify subjects connectivity matrix ...")
    print("This will take a moment, please do not interupt the process ...")
    
    # Find all subjects in phenotype file
    all_subjects = set(pheno['participant_id'])
    
    # Find subjects processed by HALFpipe
    processed_subjects = set()
    for rid, row in tqdm(pheno.iterrows()):
        file_path = glob.glob(os.path.join(derivatives_p, confounds_p.format(row['participant_id'], row['participant_id'])))
        
        if any(os.path.exists(path) for path in file_path):
            processed_subjects.add(row['participant_id'])
    
    # Find unprocessed subjects
    unprocessed_subjects = all_subjects - processed_subjects
    
    # Create report
    report = {
        "halfpipe_unprocessed": {
            "total_subjects": len(all_subjects),
            "processed_subjects": len(processed_subjects),
            "unprocessed_subjects": len(unprocessed_subjects),
            "unprocessed_subject_list": list(unprocessed_subjects)
        }
    }
    
    # Create results directory if it doesn't exist
    os.makedirs(out_p, exist_ok=True)
    
    # Save report to JSON
    report_file = os.path.join(out_p, 'cwas_subject_report.json')
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=4)
    
    # Print results
    print("\n=== Summary HALFpipe processing ===")
    print(f"Total subjects: {len(all_subjects)}")
    print(f"Processed subjects: {len(processed_subjects)}")
    print(f"Unprocessed subjects: {len(unprocessed_subjects)}")
    
    if unprocessed_subjects:
        print("\nSubjects not processed by HALFpipe:")
        for subject in sorted(unprocessed_subjects):
            print(f"- {subject}")
    else:
        print("All subjects have been processed by HALFpipe!")

    print(f"Information saved in:", report_file)
    
    return pheno[pheno['participant_id'].isin(processed_subjects)].reset_index(drop=True)