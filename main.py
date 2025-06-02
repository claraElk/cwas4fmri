import argparse
from workflow import run_pipeline

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bids_dir", required=True, help="Path to the BIDS directory")
    parser.add_argument("--output_dir", required=True, help="Path to the output directory")
    parser.add_argument("--phenotype-file", required=True, help="Path to the phenotype file")
    parser.add_argument("--atlas-file", required=True, help="Path to the atlas file")
    parser.add_argument("--analysis_level", choices=["group"])
    parser.add_argument("--participant-label", nargs="+")

    # Based on phenotype file
    parser.add_argument("--case_id", type=str, required=True, help="ID for the case group")
    parser.add_argument("--control_id", type=str, required=True, help="ID for the control group")

    # Optional arguments
    parser.add_argument("--scanner", action="store_true", default=False, help="Include scanner information")
    parser.add_argument("--sequence-col", action="store_true", default=False, help="Include sequence column in the phenotype file")
    parser.add_argument("--medication-col", action="store_true", default=False, help="Include medication column in the phenotype file")

    args = parser.parse_args()

    run_pipeline(args.bids_dir, args.output_dir, 
                 args.phenotype_file, args.atlas_file,
                 args.participant_label, args.scanner,
                 args.sequence_col, args.medication_col, 
                 args.case_id, args.control_id)
