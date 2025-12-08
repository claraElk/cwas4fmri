#!/bin/bash

# This script runs a case-wide association study (CWAS) analysis using the provided Python script.

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

for feature in "${list_features[@]}"; do

    python run_cwas.py \
        --feature_settings "$feature" \
        --pheno_filtered_path "/Users/claraelkhantour/Documents/projects/papers/halfpipe_paper2_strategies/ds000030/participants_filtered.tsv" \
        --derivatives_path "/Users/claraelkhantour/Documents/projects/papers/halfpipe_paper2_strategies/ds000030/derivatives/halfpipe" \
        --atlas_file "/Users/claraelkhantour/Documents/projects/papers/halfpipe_paper2_strategies/atlas-Schaefer2018Combined_dseg.tsv" \
        --out_p "/Users/claraelkhantour/Documents/projects/papers/halfpipe_paper2_strategies/ds000030/cwas" \
        --case_name "SCHZ" \
        --control_name "CONTROL" \
        --sequence_col "ScannerSerialNumber" \

done