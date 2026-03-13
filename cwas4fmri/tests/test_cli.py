"""
Smoke test of cwas4fmri.
"""

import pytest
import pandas as pd
import numpy as np

from cwas4fmri.run import global_parser
from pathlib import Path
from cwas4fmri.workflow import workflow

import importlib.resources as pkg_resources


@pytest.mark.smoke
def test_halfpipe(tmp_path: Path):
    data_path = Path(
        pkg_resources.files("cwas4fmri").joinpath(
            "data/test_data/dataset-ds000030_downscaled_halfpipe1.2.3dev"
        )
    )
    atlas_label = "schaefer400"
    dseg_path = data_path / "atlas" / "atlas-Schaefer2018Combined_dseg.tsv"

    bids_dir = tmp_path / data_path / "derivatives"

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    subjects = [
        f"sub-{i}"
        for i in [
            "10159",
            "10171",
            "10189",
            "10206",
            "10217",
            "10225",
            "10227",
            "10228",
            "10235",
            "10249",
        ]
    ]

    phenotypes = pd.DataFrame(
        dict(
            participant_id=subjects,
            age=np.random.uniform(18, 80, len(subjects)),
            gender=np.random.choice(["M", "F"], len(subjects)),
            diagnosis=np.random.choice(["case", "control"], len(subjects)),
            numerical_covariates=np.random.uniform(0, 10, len(subjects)),
            categorical_covariates=np.random.choice(
                ["Site1", "Site2"], len(subjects)
            ),
        )
    )
    phenotypes_path = bids_dir / "participants.tsv"
    phenotypes.to_csv(phenotypes_path, sep="\t", index=False)

    parser = global_parser()

    argv = [
        str(bids_dir),
        str(output_dir),
        "group",
        "--strategy",
        "corrMatrix1",
        "--phenotype",
        str(phenotypes_path),
        "--atlas",
        atlas_label,
        "--atlas_file",
        str(dseg_path),
        "--patient",
        "case",
        "--control",
        "control",
        "--categorical_covariates",
        "categorical_covariates",
        "--numerical_covariates",
        "numerical_covariates",
    ]

    args = parser.parse_args(argv)
    workflow(args)

    # Save results
    base_filename = (
        f"{args.patient}-{args.control}"
        f"_feature-{args.strategy}_atlas-{args.atlas}_desc-cwas"
    )

    assert (output_dir / f"{base_filename}_standardized_betas.tsv").is_file()
    assert (
        output_dir / f"{base_filename}_fdr_corrected_pvalues.tsv"
    ).is_file()
    assert (output_dir / f"{base_filename}_pvalues.tsv").is_file()
    assert (output_dir / f"{base_filename}_betas.tsv").is_file()
