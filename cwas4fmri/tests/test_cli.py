"""
Smoke test of cwas4fmri.
"""

import json
import random

import pandas as pd
import numpy as np

from cwas4fmri.run import global_parser
from pathlib import Path
from cwas4fmri.workflow import workflow


def create_fake_dataset(tmp_path: Path):
    # Create fake symmetric connectivity matrices for 10 subjects
    halfpipe_dir = tmp_path / "derivatives" / "halfpipe"
    halfpipe_dir.mkdir(parents=True, exist_ok=True)
    print(halfpipe_dir)

    n_subjects = 10
    n_rois = 5
    for i in range(n_subjects):
        subj_id = f"sub-{i:02d}"
        mat = np.random.rand(n_rois, n_rois)
        mat = (mat + mat.T) / 2  # Make it symmetric
        np.fill_diagonal(mat, 1)  # Set diagonal to 1

        # Save as TSV
        out_path = (
            halfpipe_dir
            / f"{subj_id}"
            / "func"
            / "task-rest"
            / f"{subj_id}_task-rest-ses-01_feature-test_"
            "atlas-Schaefer2018Combined_desc-correlation_matrix.tsv"
        )
        out_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(mat).to_csv(out_path, sep="\t", index=False, header=False)

        # Create corresponding JSON with FDMean and FDMax
        json_path = (
            halfpipe_dir
            / f"{subj_id}"
            / "func"
            / "task-rest"
            / f"{subj_id}_task-rest-ses-01_feature-test_"
            "atlas-Schaefer2018Combined_timeseries.json"
        )
        json_data = {
            "FDMean": random.randrange(0, 1),
            "FDMax": random.randrange(0, 4),
        }
        with open(json_path, "w") as f:
            json.dump(json_data, f)

    # Create a participants.tsv file
    participants_path = tmp_path / "participants.tsv"
    participant_data = {
        "participant_id": [f"sub-{i:02d}" for i in range(n_subjects)],
        "diagnosis": [
            "SCHZ" if i < n_subjects // 2 else "CONTROL"
            for i in range(n_subjects)
        ],
        "age": np.random.randint(20, 60, size=n_subjects),
        "gender": ["M" if i % 2 == 0 else "F" for i in range(n_subjects)],
    }
    pd.DataFrame(participant_data).to_csv(
        participants_path, sep="\t", index=False
    )

    # Create a fake atlas dseg file
    atlas_dir = tmp_path / "atlases"
    atlas_dir.mkdir(parents=True, exist_ok=True)
    dseg_path = atlas_dir / "atlas-Schaefer2018Combined_dseg.tsv"
    dseg_data = np.arange(1, n_rois + 1).reshape(-1, 1)
    pd.DataFrame(dseg_data).to_csv(
        dseg_path, sep="\t", index=True, header=False
    )


def test_cli(tmp_path: Path):

    create_fake_dataset(tmp_path)

    atlas_label = "Schaefer2018Combined"
    dseg_path = tmp_path / "atlases" / "atlas-Schaefer2018Combined_dseg.tsv"

    bids_dir = tmp_path / "derivatives"

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    phenotypes_path = tmp_path / "participants.tsv"

    parser = global_parser()

    argv = [
        str(bids_dir),
        str(output_dir),
        "group",
        "--strategy",
        "test",
        "--phenotype",
        str(phenotypes_path),
        "--atlas",
        atlas_label,
        "--atlas_file",
        str(dseg_path),
        "--patient",
        "SCHZ",
        "--control",
        "CONTROL",
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
