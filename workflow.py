from tqdm import tqdm

from utils.phenotype import load_phenotype
from utils.files import verify_files_and_directories, verify_atals_files, find_halfpipe_output
from utils.subject import find_valid_subjects
from utils.reject_fd_qc import filter_by_qc, filter_by_fd
from utils.features import extract_feature_settings
from utils.connectome import process_connectivity_matrix
from utils.stats import define_regressors, glm_wrap_cc, summarize_glm, save_glm

def run_pipeline(working_directory, output_dir, atlas_file, pheno_p,
                 scanner, sequence_col, medication_col, 
                 case_name, control_name,) : 

    dict_halfpipe = find_halfpipe_output(working_directory)

    # 1. Verify the phenotype file
    verified_df = load_phenotype(pheno_p, 
                                 diagnosis_col="diagnosis", # Hardcoded required columns
                                 subject_col="participant_id", 
                                 age_col="age", 
                                 sex_col="sex", 
                                 scanner_col=scanner,
                                 sequence_col=sequence_col, medication_col=medication_col,
                                 case_name=case_name, control_name=control_name,
                                 subject_file_path=dict_halfpipe['subject_p']
                                 )

    # 2. Verify atlas location and format: only accept tsv file
    labels, conn_mask, roi_labels = verify_atals_files(atlas_file)

    print(dict_halfpipe['json_exclude_qc_path'])
    verify_files_and_directories(working_directory, 
                                dict_halfpipe['json_exclude_qc_path'], dict_halfpipe['reports_dir'])

    # 3. Find number of subjects processed by HALFpipe, in case some subjects failed
    pheno_filtered = find_valid_subjects(verified_df, derivatives_p=dict_halfpipe['derivatives_p'],
                                        confounds_p=dict_halfpipe['confounds_p'],
                                        out_p=output_dir
                                        )
    
    # 4. Reject based on bad QC & FD threshold
    #pheno_filtered_qc = filter_by_qc(dict_halfpipe['json_exclude_qc_path'], pheno_filtered, out_p=output_dir)

    # 5. Reject subject based on FD>0.5
    pheno_filtered_qc_fd = filter_by_fd(
                            pheno_filtered_qc=pheno_filtered, # change pheno_filtered with pheno_filtered_qc if you have run step 4.
                            derivatives_p=dict_halfpipe['derivatives_p'],
                            confounds_p=dict_halfpipe['confounds_p'],
                            out_p=output_dir
                            )

    # Extract feature settings and validate atlases
    feature_settings, common_atlas = extract_feature_settings(dict_halfpipe['json_spec_path'])
        
    # Define regressors
    regressors = define_regressors(scanner, sequence_col, medication_col)

    # Process each feature
    for feature in tqdm(feature_settings):
        print(f"\nProcessing feature: {feature}")
        
        # Process connectivity matrix
        conn_stack, final_df = process_connectivity_matrix(
            pheno_filtered_fd=pheno_filtered_qc_fd,
            connectome_t=dict_halfpipe['connectome_t'],
            feature=feature,
            atlas=common_atlas,
            connectome_p=dict_halfpipe['connectome_p'],
            conn_mask=conn_mask
        )
        
        # Perform CWAS analysis
        glm_con = glm_wrap_cc(conn_stack, final_df, 
                                group='diagnosis', case=1, control=0, 
                                regressors=regressors, report=True)
        
        # Get results
        table_con, table_stand_beta_con, table_qval_con = summarize_glm(
            glm_con, conn_mask, roi_labels
        )


        save_glm(output_dir, table_con, 
                table_stand_beta_con, table_qval_con, 
                conn_mask, roi_labels,
                case_name, control_name, 
                feature, common_atlas)

