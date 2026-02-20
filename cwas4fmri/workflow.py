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
        args (argparse.Namespace): Arguments parsed from the global parser

    Returns:
        None
    """
    set_verbosity(args.verbosity)

    logger.info(
        f"Run CWAS analysis for :"
        f"\n-strategy: {args.strategy}"
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
    results = summarize_glm(glm_con, mask_2d, roi_labels)

    # Save results
    base_filename = (
        f"{args.patient}-{args.control}"
        f"_feature-{args.strategy}_atlas-{args.atlas}_desc-cwas"
    )

    results["out_table"].to_csv(output / f"{base_filename}.tsv", sep="\t")
    results["stand_beta_table"].to_csv(
        output / f"{base_filename}_standardized_betas.tsv", sep="\t"
    )
    results["qval_table"].to_csv(
        output / f"{base_filename}_fdr_corrected_pvalues.tsv", sep="\t"
    )
    results["pval_table"].to_csv(
        output / f"{base_filename}_pvalues.tsv", sep="\t"
    )
    results["beta_table"].to_csv(
        output / f"{base_filename}_betas.tsv", sep="\t"
    )

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
