from cwas4fmri.phenotype import load_phenotype
from cwas4fmri.subject import find_valid_subjects
from cwas4fmri.reject_fd_qc import filter_by_fd 
from cwas4fmri.connectome import process_connectivity_matrix
from cwas4fmri.plots import plot_interactive_matrix
from cwas4fmri.stats import *
from cwas4fmri.files import *

def run_pipeline(bids_dir, output_dir, pheno_p, atlas_file, atlas, group,
                scanner, sequence, medication, case_name, control_name, 
                session, task, run, feature): 

    bids_validation(bids_dir=bids_dir)

    create_output_directory(output_dir)

    dict_halfpipe = find_bids_output(bids_dir)

    # Verify the phenotype file
    df = load_phenotype(pheno_p,
                        diagnosis_col=group,
                        subject_col="participant_id",
                        age_col="age",
                        sex_col="sex",
                        scanner_col=scanner,
                        sequence=sequence, medication=medication,
                        case_name=case_name, control_name=control_name,
                        )
  
    # Verify atlas location and format
    conn_mask, roi_labels = verify_atlas_files(atlas_file)

    # Find number of subjects
    df_filtered = find_valid_subjects(
        bids_dir=bids_dir,
        pheno=df,
        session=session,
        connectome_t=dict_halfpipe['connectome_t'],
        run=run,
        task=task,
        atlas=atlas,
        feature=feature,
        out_p=output_dir
        )

    # Reject subject based on FD>0.5
    pheno_filtered_qc_fd = filter_by_fd(
        pheno_filtered_qc=df_filtered,
        derivatives_p=bids_dir,
        confounds_json=dict_halfpipe['confounds_json'],
        out_p=output_dir,
        session=session,
        task=task,
        run=run,
        feature=feature
        )

    # Define regressors
    regressors = define_regressors(scanner, sequence, medication)

    # Process connectivity matrix
    conn_stack, final_df = process_connectivity_matrix(
        pheno_filtered_fd=pheno_filtered_qc_fd,
        connectome_t=dict_halfpipe['connectome_t'],
        feature=feature,
        atlas=atlas,
        bids_dir=bids_dir,
        conn_mask=conn_mask,
        session=session, 
        task=task,
        run=run
        )

    # Perform CWAS analysis
    glm_con = glm_wrap_cc(output_dir, conn_stack, final_df,
                            group=group, case=1, control=0, 
                            regressors=regressors, report=True)

    # Get results
    table_con, table_stand_beta, table_qval = summarize_glm(
        glm_con, 
        conn_mask, 
        roi_labels
        )

    save_glm(
        out_p=output_dir, 
        table_con=table_con,
        table_stand_beta_con=table_stand_beta,
        table_qval_con=table_qval,
        conn_mask=conn_mask,
        roi_labels=roi_labels,
        case_name=case_name,
        control_name=control_name,
        feature=feature,
        atlas=atlas)
    
    plot_interactive_matrix(output_path=output_dir,
                            beta_matrix=table_stand_beta, 
                            pvalues_matrix=table_qval, 
                            labels=roi_labels
                            )