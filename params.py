import os
from pathlib import Path

# Configuration
working_directory = "path/to/working_dir" # Path to your working directory
atlas_file = "path/to/atlas.tsv" # Path to acces atlas
pheno_p = 'path/to/pheno.tsv'  # Path to phenotype file
out_p = 'path/to/output' # Path to results, where results will be saved, the file will be automatically created


#################################################################################
# DO NOT MODIFY THE SECTION BELOW
# all your files will be automatically detected.
# However, if you moved filed from the original working directory, you can indicate the paths below 

reports_dir = os.path.join(working_directory, "reports") # Path to report folder
subject_p = os.path.join(working_directory, "subject-list.txt") # Subject list

derivatives_p = os.path.join(working_directory, "derivatives") # Path to derivative folder
connectome_p = os.path.join(derivatives_p, 'halfpipe') # Path to halfpipe derivative folder

# Create output folder
out_p = Path(out_p)
print("Creating output folder if does not exist ...")
print(out_p)
out_p.mkdir(parents=True, exist_ok=True) # Create folder if doesn't exist

json_spec_path = os.path.join(working_directory, "spec.json")  # Path to spec.json file
json_exclude_qc_path = os.path.join(reports_dir, 'exclude.json') # Path to exclude.json file

def check_path():
    connectome_t = os.path.join('{}', '*func', '*', '{}_*_feature-{}_atlas-{}_desc-correlation_matrix.tsv')
    confounds_p = os.path.join('fmriprep', '{}', '*func', '*{}_*_desc-confounds_timeseries.tsv')
    return connectome_t, confounds_p