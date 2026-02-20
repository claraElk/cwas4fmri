# 🧠 CWAS for resting-state fMRI

This workflow is compliant with atlas-based connectomes from HALFpipe v1.3.0 or higher

# How to run cwas4fmri
1. Create an environment
2. Install the package:
```bash
pip install .
```
3. Then run the pipeline
```bash
cwas4fmri bids_dir \
          output_dir \
          group \
          --strategy \
          --phenotype \
          --atlas \
          --atlas_file \
          --patient \
          --control \
          --categorical_covariates \
          --numerical_covariates \
          --verbosity \
          --debug
```

For more details, run `cwas4fmri -h`

# HALFpipe folder and participants.tsv file
> **Note:** The pipeline processes only subjects present in **both** the phenotype file and the HALFpipe derivatives folder.
> Ensure both are clean before running:
> - Subjects absent from the phenotype file will be skipped, even if their connectomes exist.
> - Subjects with a missing connectome will be excluded from the analysis.
> *The phenotype file is the reference:* only subjects listed there will be processed.

- Your participants.tsv file need to have the following columns:

| participant_id | diagnosis | gender | age |
| ---------------|-----------|--------|-----|

# Environment
The environment will be automatically created and deleted with the job.
