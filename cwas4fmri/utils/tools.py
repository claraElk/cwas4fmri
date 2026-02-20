# Credit original code: Sebastian Urchs

import json
from pathlib import Path
from tqdm import tqdm
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests as stm

from .stats import conn2mat
from ..logger import logger


def reject_fd(
    json_files: list, phenotype: pd.DataFrame, subj: str
) -> pd.DataFrame:
    """
    Extract the mean FD in HALFpipe json file and average through runs

    Args:
        json_files (list): list of path to 1 subject json file(s)

    Returns:
        np.float() : averaged mean_fd across runs
    """

    exclude_subject = []

    for jf in json_files:
        with open(jf, "r") as f:
            metadata = json.load(f)
        if metadata["FDMax"] > 3.0:
            exclude_subject.append(subj)
        if metadata["FDMean"] > 0.5:
            exclude_subject.append(subj)

    return phenotype[
        ~phenotype.participant_id.isin(exclude_subject)
    ].reset_index(drop=True)


def summarize_glm(
    glm_table: pd.DataFrame, mask_2d: np.ndarray, labels: list
) -> pd.DataFrame:
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
    beta_table = pd.DataFrame(
        conn2mat(out_table.betas.values, mask_2d), index=labels, columns=labels
    )

    stand_beta_table = pd.DataFrame(
        conn2mat(out_table.stand_betas.values, mask_2d),
        index=labels,
        columns=labels,
    )

    pval_table = pd.DataFrame(
        conn2mat(out_table.pvals.values, mask_2d), index=labels, columns=labels
    )

    qval_table = pd.DataFrame(
        conn2mat(out_table.qval.values, mask_2d), index=labels, columns=labels
    )

    return out_table, stand_beta_table, qval_table, pval_table, beta_table


def average_runs(corr_files: list) -> np.ndarray:
    """
    Average connectomes through runs
    Args:
        corr_files (list): list of path of HALFpipe correlation matrices

    Returns:
        avg_mat (np): Averaged connectivity matrix of 1 subject (n roi x n roi)
    """

    # Load all matrices for this subject
    matrices = [
        pd.read_csv(cf, sep="\t", header=None, dtype=np.float32).values
        for cf in corr_files
    ]

    # Average across runs using nanmean
    avg_mat = np.nanmean(matrices, axis=0)

    return avg_mat


def fd_mean_extraction(json_files: list, subj: str) -> float:
    """
    Extract the mean FD in HALFpipe json file and average through runs

    Args:
        json_files (list): list of path to 1 subject json file(s)

    Returns:
        np.float() : averaged mean_fd across runs
    """

    fd_values = []

    for jf in json_files:
        with open(jf, "r") as f:
            metadata = json.load(f)
        if "FDMean" in metadata:
            fd_values.append(metadata["FDMean"])
        else:
            logger.warning(f"Warning: could not read JSON for {subj}")

    return np.mean(fd_values) if fd_values else np.nan


def process_connectivity_matrix(
    phenotype: pd.DataFrame, feature: str, derivatives_path: str
) -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    """
    Process connectivity matrices and return CWAS-ready flattened matrices.
    - Averages multiple runs per subject
    - Keeps NaNs
    - Uses upper triangle without diagonal

    Args:
        phenotype (pd.DataFrame) : dataframe with only good QC subjects
            - Columns age, gender and diagnosis are mandatory
        feature (str): pipeline name
        derivative_path (str): path to HALFpipe derivatives
    """

    logger.info("Processing connectivity matrices ...")

    flattened_matrices = []
    valid_subject_indices = []
    fdmean_values = []

    # Reject subjects based on FD
    json_files = list(
        Path(derivatives_path).rglob(
            f"sub-*_feature-{feature}_*_timeseries.json"
        )
    )
    phenotype = reject_fd(
        json_files=json_files, phenotype=phenotype, subj=None
    )

    for idx, row in tqdm(phenotype.iterrows(), total=len(phenotype)):
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
        fdmean_values.append(fd_mean_extraction(json_files, subj))

    # --- Build final conn_stack ---
    conn_stack = np.vstack(flattened_matrices)  # shape: (n_subjects, n_edges)

    # TODO: WARNING, this is tricky, should be changed!!!
    phenotype = phenotype.loc[valid_subject_indices].reset_index(drop=True)
    phenotype["mean_fd"] = fdmean_values

    logger.info(f"Processed {len(phenotype)} subjects")
    logger.info(
        f"Connectivity stack shape (n_subjects x n_edges): {conn_stack.shape}"
    )

    return conn_stack, phenotype, mask_2d
