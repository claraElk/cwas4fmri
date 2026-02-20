#!/bin/bash
#SBATCH --job-name=cwas4fmri
#SBATCH --output=logs/cwas4fmri_%A_%a.log
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=0-17   # <-- number of features - 1

list_features=(
  'Baseline' 
  'Wang2023Simple' 
  'Wang2023SimpleGSR'
  'Wang2023ScrubbingGSR'
  'Wang2023Scrubbing'
   'Wang2023aCompCorGSR'
   'Wang2023aCompCor'
   'motionParametersScrubbingGSR'
   'motionParametersScrubbing'
   'motionParametersGSR'
   'motionParameters'
   'ICAAROMAScrubbingGSR'
   'ICAAROMAScrubbing'
   'ICAAROMAGSR'
   'ICAAROMACCompcor'
   'ICAAROMA'
   'Everything'
   'cCompCor'
)

feature=${list_features[$SLURM_ARRAY_TASK_ID]}

echo "Running feature: $feature"

cd path/to/cwas4fmri

module load python/3.11
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index nilearn pandas "numpy<2" scipy statsmodels argparse scikit-learn tqdm

python run_cwas.py \
    --feature_settings "$feature" \
    --pheno_filtered_path path/to/phenotype.tsv \
    --derivatives_path path/to/derivatives \
    --atlas_file path/to/atlas.tsv \
    --out_p path/to/cwas \
    --case_name "case" \
    --control_name "Control" \
    --sequence_col "sequence"
