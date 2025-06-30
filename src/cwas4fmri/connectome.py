import os
from tqdm import tqdm
import numpy as np
import pandas as pd

def conn2mat(conn, mask):
    conn_mat = np.zeros(shape=mask.shape)
    conn_mat[mask] = conn
    conn_mat += conn_mat.T
    conn_mat[np.eye(mask.shape[0]).astype(bool)] /= 2

    return conn_mat

def process_connectivity_matrix(pheno_filtered_fd, connectome_t, feature, atlas, bids_dir, conn_mask, session, task, run):
    """
    Process connectivity matrices for valid subjects
    """
    print("\n⏳ Process connectivity matrices for valid subjects ...")

       # Collect valid connectome paths
    valid_subject_paths = []
    valid_subject_indices = []
    
    for index, row in tqdm(pheno_filtered_fd.iterrows()):
        connectome_file = os.path.join(
            bids_dir,
            connectome_t.format(
                row['participant_id'], session, 
                row['participant_id'], session, 
                task, run, atlas, feature
            )
        )
        if os.path.exists(connectome_file):  # Check if file exists
            valid_subject_paths.append(connectome_file)
            valid_subject_indices.append(index)  # Store index of valid subjects

    # Stack connectome data
    conn_stack = np.array([pd.read_csv(p, sep='\t').values[conn_mask] for p in valid_subject_paths])
    
    print(f"\n📌 Processing statistics for feature {feature} with atlas {atlas}:")
    
    return conn_stack, pheno_filtered_fd