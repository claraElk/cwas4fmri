import os
import glob
from tqdm import tqdm
import numpy as np
import pandas as pd

def conn2mat(conn, mask):
    conn_mat = np.zeros(shape=mask.shape)
    conn_mat[mask] = conn
    conn_mat += conn_mat.T
    conn_mat[np.eye(mask.shape[0]).astype(bool)] /= 2

    return conn_mat

def process_connectivity_matrix(pheno_filtered_fd, connectome_t, feature, atlas, connectome_p, conn_mask):
    """
    Process connectivity matrices for valid subjects
    
    Args:
        pheno_filtered_fd (pd.DataFrame): Filtered phenotype data with FD < 0.5
        feature (str): Feature name
        atlas (str): Atlas name
        connectome_p (str): Path to connectome directory
        conn_mask (array): Mask for connectome data
    
    Returns:
        tuple: (conn_stack, pheno_filtered_fd)
    """
    print("Process connectivity matrices for valid subjects ...")

    # Collect valid connectome paths
    valid_subject_paths = []
    valid_subject_index = []
    valid_subject_ids = []
    
    for index, subject in tqdm(pheno_filtered_fd.iterrows()):
        connectome_file = glob.glob(os.path.join(connectome_p, connectome_t.format(
            subject['participant_id'],subject['participant_id'], feature, atlas
        )), recursive=True)
        
        for file in connectome_file : 
            if os.path.exists(file):  # Check if file exists
                valid_subject_paths.append(file)
                valid_subject_index.append(index)  # Store index of valid subjects
                valid_subject_ids.append(subject['participant_id'])  # Track subject ID
     
    # Stack connectome data
    conn_stack = np.array([pd.read_csv(p, sep='\t', header=None).values[conn_mask] for p in valid_subject_paths])
    
    conn_df = pd.DataFrame(conn_stack)
    conn_df['subject_id'] = valid_subject_ids
    
    # Averaged per subject
    conn_avg = conn_df.groupby('subject_id').mean()

    # Align conn with pheno
    # Ensure pheno has subject_id as index
    pheno_aligned = pheno_filtered_fd.set_index('participant_id').loc[conn_avg.index]

    # Get statistics
    stats = {
        'subjects': len(pheno_aligned),
        'valid_connectome': len(valid_subject_paths)
    }
    
    print(f"\nProcessing statistics for feature {feature} with atlas {atlas}:")
    
    return conn_avg.values, pheno_aligned