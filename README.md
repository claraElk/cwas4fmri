# 🧠 CWAS for resting-state fMRI

This workflow is compliant with atlas-based connectomes from HALFpipe v1.3.0 or higher

# How to run cwas4fmri
1. `git clone https://github.com/claraElk/cwas4fmri.git`
2. Create a virtual environment
3. Install the package:
```bash
pip install .
```
4. Then run the pipeline
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

# How to run cwas4fmri with docker (not fully tested yet)
1. `git clone https://github.com/claraElk/cwas4fmri.git`
2. `docker build -t cwas4fmri:latest`
3. Run using docker
```
docker run -v /Users/username/:/Users/username/ cwas4fmri:latest bids_dir \
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
*Note: To access your files more easily, you can bind in the path where the working directory is located with the `-v` flag*
