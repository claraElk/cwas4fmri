import argparse
from cwas4fmri.workflow import run_pipeline

def parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bids_dir", required=True, help="Path to the BIDS directory")
    parser.add_argument("--output_dir", required=True, help="Path to the output directory")
    parser.add_argument("--phenotype_file", required=True, help="Path to the phenotype file")
    
    parser.add_argument("--atlas", required=True, help="Atlas to use for the analysis")
    parser.add_argument("--atlas_file", required=True, help="Path to the atlas file")

    parser.add_argument("--analysis_level", choices=["group"])
    parser.add_argument("--participant-label", nargs="+")
    
    parser.add_argument("--session", type=str, required=True, help="Session label for the analysis")
    parser.add_argument("--task", type=str, required=True, help="Task label for the analysis")
    parser.add_argument("--run", type=str, required=True, help="Run label for the analysis")
    parser.add_argument("--feature", type=str, required=True, help="Feature label for the analysis")

    # Based on phenotype file
    parser.add_argument("--group", type=str, required=False, default="diagnosis", help="Column name for the diagnosis in the phenotype file")
    parser.add_argument("--case_id", type=str, required=True, help="ID for the case group")
    parser.add_argument("--control_id", type=str, required=True, help="ID for the control group")

    # Optional arguments
    parser.add_argument("--scanner", type=bool, default=False, help="Include scanner information")
    parser.add_argument("--sequence", type=bool, default=False, help="Include sequence column in the phenotype file")
    parser.add_argument("--medication", type=bool, default=False, help="Include medication column in the phenotype file")

    args = parser.parse_args()

    return args

def main():
    print("\n🚀 Welcome to CWAS-rsfmri! \n")
    args = parsers()

    run_pipeline(bids_dir=args.bids_dir, 
                 output_dir=args.output_dir, 
                 pheno_p=args.phenotype_file, 
                 atlas_file=args.atlas_file,
                 atlas=args.atlas, 
                 group=args.group,
                 scanner=args.scanner,
                 sequence=args.sequence, 
                 medication=args.medication, 
                 case_name=args.case_id, 
                 control_name=args.control_id, 
                 session=args.session,
                 task=args.task, 
                 run=args.run, 
                 feature=args.feature)
    
    print("\n🎉 Pipeline finished! \n")
