# Credit original code: Sebastian Urchs

import json
from pathlib import Path
from tqdm import tqdm
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests as stm

from .stats import conn2mat
from ..logger import logger


def filter_and_extract_fd(json_files: list, subj: str) -> tuple[float | None]:
    """
    For a single subject, check FD criteria across all their runs and
    extract their mean FD value.
    Args:
        json_files (list): List of paths to JSON file(s) for this subject.
        subj (str): Subject identifier.
    Returns:
        tuple:
            - bool: True if subject passes FD criteria, False otherwise.
            - float | None: Mean FD averaged across runs, or None if excluded.
    """
    fdmean_values = []

    for jf in json_files:
        with open(jf, "r") as f:
            metadata = json.load(f)

        if not metadata:
            logger.warning(
                f"Could not read JSON for {subj} at {jf}, excluding"
            )
            return None

        if metadata["FDMax"] > 3.0:
            logger.warning(f"Excluding {subj} based on FD max criteria")
            return None

        if metadata["FDMean"] > 0.5:
            logger.warning(f"Excluding {subj} based on FD mean criteria")
            return None

        fdmean_values.append(metadata["FDMean"])

    return float(np.mean(fdmean_values))


def summarize_glm(
    glm_table: pd.DataFrame, mask_2d: np.ndarray, labels: list
) -> pd.DataFrame:
    """
    Summarize GLM results:
    - Converts flattened upper-triangle edges back to full square matrix
    - Computes FDR q-values
    Args:
        glm_table (pd.DataFrame): DataFrame with results for each edge
        mask_2d (np.ndarray): Indicating upper-triangle positions
        labels (list): List of ROI labels

    Returns:
        dict: Dictionary containing:
            - out_table: Original glm_table with added 'qval' column
            - beta_table: DataFrame of full square matrix of betas
            - stand_beta_table: DataFrame of standardized betas
            - pval_table: DataFrame of full square matrix of p-values
            - qval_table: DataFrame of full square matrix of q-values
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

    results = {
        "out_table": out_table,
        "stand_beta_table": stand_beta_table,
        "qval_table": qval_table,
        "pval_table": pval_table,
        "beta_table": beta_table,
    }
    return results


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
    Returns:
        conn_stack (np.ndarray): 2D array of shape (n_subjects, n_edges)
        phenotype (pd.DataFrame): Updated phenotype dataframe and added mean FD
        mask_2d (np.ndarray): Positions for reconstructing full matrices
    """
    logger.info("Processing connectivity matrices ...")

    records = []

    for _, row in tqdm(phenotype.iterrows(), total=len(phenotype)):
        subj = str(row["participant_id"])
        if not subj.startswith("sub-"):
            subj = f"sub-{subj}"

        # --- Check and filter subject based on FD criteria ---
        json_pat = f"{subj}_*feature-{feature}_*_timeseries.json"
        json_files = list(Path(derivatives_path).rglob(json_pat))
        if not json_files:
            logger.warning(
                f"No JSON file found for {subj}, excluding from analysis"
            )
            continue

        fdmean = filter_and_extract_fd(json_files, subj)
        if fdmean is None:
            continue

        # --- Check connectivity matrix files exist ---
        corr_pat = f"{subj}_*feature-{feature}_*desc-correlation_matrix.tsv"
        corr_files = list(Path(derivatives_path).rglob(corr_pat))
        if not corr_files:
            logger.warning(
                f"No correlation matrix found for {subj},"
                "excluding from analysis"
            )
            continue

        # --- Average runs ---
        avg_mat = average_runs(corr_files)
        n_rois = avg_mat.shape[0]

        # TODO: replace with nilearn mat_to_vec function
        mask_2d = np.triu(np.ones((n_rois, n_rois), dtype=bool), k=1)

        # All checks passed: store result explicitly based on participant_id
        records.append(
            {
                "participant_id": row[
                    "participant_id"
                ],  # use original id from phenotype
                "flattened": avg_mat[mask_2d],
                "mean_fd": fdmean,
            }
        )

    if not records:
        raise ValueError("No valid subjects remaining after processing.")

    conn_stack = np.vstack(
        [r["flattened"] for r in records]
    )  # shape: (n_subjects, n_edges)

    # Merge FDmean by participant_id
    fd_df = pd.DataFrame(
        {
            "participant_id": [r["participant_id"] for r in records],
            "mean_fd": [r["mean_fd"] for r in records],
        }
    )
    valid_ids = fd_df["participant_id"].tolist()
    phenotype = phenotype[
        phenotype["participant_id"].isin(valid_ids)
    ].reset_index(drop=True)
    phenotype = phenotype.merge(fd_df, on="participant_id", how="left")

    logger.info(f"Processed {len(phenotype)} subjects")
    logger.info(
        f"Connectivity stack shape (n_subjects x n_edges): {conn_stack.shape}"
    )

    return conn_stack, phenotype, mask_2d
