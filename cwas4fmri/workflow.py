from .utils.tools import summarize_glm, process_connectivity_matrix
from .utils.stats import glm_wrap_cc
from pathlib import Path
import pandas as pd
import argparse
from .logger import logger, set_verbosity


def workflow(args):
    """
    Run CWAS analysis for all features defined in the JSON specification

    Args:
        feature (str): Pipeline name (as indicated in HALFpipe filename)
        phenotype (pd.DataFrame): Filtered phenotype data (only good subjects)
        derivatives_p (str): Path to derivatives directory
        roi_labels (list): List of ROI labels
        output (str): Output directory path
        patient (str): case variable as indicated in phenotype file
        control (str): control group variable as indicated in phenotype file
        cat_covariates (str): list of categorical covariates
        num_covariates (str): list of numerical covariates
        atlas (str): name of the atlas, as indicated in HALFpipe filename.
            - Default: schaeferCombined

    Returns:
        None
    """
    set_verbosity(args.verbosity)

    logger.info(
        f"Run CWAS analysis for :"
        f"\n-strategy {args.strategy}"
        f"\n-atlas: {args.atlas}"
        f"\n-patient group: {args.patient}"
        f"\n-control group: {args.control}"
    )

    # Create output directory if it doesn't exist
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)

    labels = pd.read_csv(args.atlas_file, sep="\t", header=None)
    roi_labels = labels[1].to_list()

    # Define regressors
    list_regressor = ["age", "mean_fd", "C(gender)"]
    if args.categorical_covariates:
        for cov in args.categorical_covariates.split(","):
            list_regressor.append(f"C({cov.strip()})")
    if args.numerical_covariates:
        for cov in args.numerical_covariates.split(","):
            list_regressor.append(cov.strip())

    regressors = " + ".join(list_regressor)
    logger.info(f"Regressors used in the model: {regressors}")

    pheno = load_data_frame(args)

    # Process connectivity matrix
    conn_stack, phenotype_cwas, mask_2d = process_connectivity_matrix(
        phenotype=pheno,
        feature=args.strategy,
        derivatives_path=args.bids_dir,
    )

    # Perform CWAS analysis
    glm_con = glm_wrap_cc(
        conn_stack,
        phenotype_cwas,
        group="diagnosis",
        case=args.patient,
        control=args.control,
        regressors=regressors,
        report=True,
    )

    # Get results
    # TODO: change this with table list or dict
    (
        table_con,
        table_stand_beta_con,
        table_qval_con,
        table_pval_con,
        table_beta_con,
    ) = summarize_glm(glm_con, mask_2d, roi_labels)

    # Save results
    base_filename = (
        f"{args.patient}-{args.control}"
        f"_feature-{args.strategy}_atlas-{args.atlas}_desc-cwas"
    )

    table_con.to_csv(output / f"{base_filename}.tsv", sep="\t")
    table_stand_beta_con.to_csv(
        output / f"{base_filename}_standardized_betas.tsv", sep="\t"
    )
    table_qval_con.to_csv(
        output / f"{base_filename}_fdr_corrected_pvalues.tsv", sep="\t"
    )
    table_pval_con.to_csv(output / f"{base_filename}_pvalues.tsv", sep="\t")
    table_beta_con.to_csv(output / f"{base_filename}_betas.tsv", sep="\t")

    logger.info(
        f"Completed processing for the following strategy: {args.strategy}"
    )
    logger.info(f"Analysis complete. Results saved to: {output}")


def load_data_frame(args: argparse.Namespace) -> pd.DataFrame:
    """Load a phenotype TSV"""
    data_frame = pd.read_csv(
        args.phenotype,
        sep="\t",
        dtype={"participant_id": str},
    )
    if "gender" not in data_frame.columns:
        raise ValueError('Phenotypes file is missing the "gender" column')
    if "age" not in data_frame.columns:
        raise ValueError('Phenotypes file is missing the "age" column')
    if "diagnosis" not in data_frame.columns:
        raise ValueError('Phenotypes file is missing the "diagnosis" column')
    return data_frame
