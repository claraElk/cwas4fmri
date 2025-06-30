import os
import json
import numpy as np
import patsy as pat
import pandas as pd
import statsmodels.api as sm
from sklearn import preprocessing as skp
from statsmodels.sandbox.stats.multicomp import multipletests as stm

from .connectome import conn2mat
from .subject import find_subset
from .files import report_file

def define_regressors(scanner, sequence_col, medication_col):
    """
    Define the list of regressors based on the provided columns.
    """
    list_regressor = ['age', 'C(sex)', 'mean_fd'] # Mandatory regressors
    if scanner : 
        list_regressor.append('C(scanner)')
    if sequence_col :
        list_regressor.append('C(sequence)')
    if medication_col  :
        list_regressor.append('C(medication)')
    
    regressors = ' + '.join(list_regressor)
    print(f'\n📌 regressors used in the model: {regressors}')
    return regressors


def save_glm(out_p, table_con, table_stand_beta_con, 
             table_qval_con, conn_mask, roi_labels,
             case_name, control_name, feature, atlas):
    
    out_table = table_con.copy()

    (fdr_pass, qval, _, _) = stm(table_con.pvals, alpha=0.05, method='fdr_bh')
    out_table['qval'] = qval
    
    # Save results
    base_filename = f'cwas_{case_name}_{control_name}_rsfmri_{feature}_{atlas}'
    table_con.to_csv(os.path.join(out_p, f'{base_filename}.tsv'), sep='\t')
    table_stand_beta_con.to_csv(os.path.join(out_p, f'{base_filename}_standardized_betas.tsv'), sep='\t')
    table_qval_con.to_csv(os.path.join(out_p, f'{base_filename}_fdr_corrected_pvalues.tsv'), sep='\t')
    
    print(f"\n✅ Completed processing for feature: {feature}")
    print(f"✅ Results saved to: {os.path.join(out_p, f'{base_filename}.tsv')}")
    

def summarize_glm(glm_table, conn_mask, roi_labels):
    out_table = glm_table.copy()
    (fdr_pass, qval, _, _) = stm(glm_table.pvals, alpha=0.05, method='fdr_bh')
    out_table['qval'] = qval

    # Return to matrix form
    stand_beta_table = pd.DataFrame(conn2mat(out_table.stand_betas.values, conn_mask) , index=roi_labels, columns=roi_labels)
    qval_table = pd.DataFrame(conn2mat(out_table.pvals.values, conn_mask), index=roi_labels, columns=roi_labels)

    return out_table, stand_beta_table, qval_table


def standardize(data, mask):
    scaler = skp.StandardScaler(with_mean=False, with_std=True)
    scaler.fit(data[mask, :])
    standardized_data = scaler.transform(data)
    return standardized_data


def find_contrast(design_matrix, contrast):
    # Find the contrast column
    contrast_columns = [(col_id, col) for col_id, col in enumerate(design_matrix.columns) if f'{contrast}' in col]
    if not len(contrast_columns) == 1:
        raise Exception(f'❌ There is no single factor that matches {contrast}: {(list(design_matrix.columns))}')
    return contrast_columns


def glm(data, design_matrix, contrast):
    contrast_id, _ = find_contrast(design_matrix, contrast)[0]
    n_data = data.shape[1]

    # Conduct the GLM
    betas = np.zeros(shape=n_data)
    pvals = np.zeros(shape=n_data)
    for conn_id in range(n_data):
        model = sm.OLS(data[:, conn_id], design_matrix)
        results = model.fit()
        betas[conn_id] = results.params.iloc[contrast_id]
        pvals[conn_id] = results.pvalues.iloc[contrast_id]

    return betas, pvals

def glm_wrap_cc(out_p, conn, pheno, group, case, control, regressors='', report=False):
    # Make sure pheno and conn have the same number of cases
    if not conn.shape[0] == pheno.shape[0]:
        print(f'❌ Connectivity matrix ({conn.shape[0]}) and phenotype file ({pheno.shape[0]}) must be same number of cases')

    # Define the subset of the sample
    sub_mask, case_masks = find_subset(pheno, group, [case, control])
    sub_conn = conn[sub_mask, :]
    sub_pheno = pheno.loc[sub_mask]
    n_sub = np.sum(sub_mask)
    n_case = np.sum(case_masks[case])
    n_control = np.sum(case_masks[control])
    n_data = sub_conn.shape[1]
    
    if report:
        summary_data = {
              f'Selected sample based on group variable': f'{group}',
              f'cases {case}': f'n={n_case}',
              f'controls {control}': f'n={n_control}',
              f'original sample': f'n={pheno.shape[0]}',
              f'new sample': f'n={n_sub}',
              f'data points available': f'{n_data}',
              f'standardized estimators are based on {group}': f'{control}'
              }
        report_file(out_p, summary_data)

        print(f'\n⏳ Performing CWAS. This might take few minutes.\n')
        for key, value in summary_data.items():
            if isinstance(value, list):
                print(f"{key}: {len(value)}")  # Or print all items if preferred
            else:
                print(f"{key}: {value}")
        
    # Standardize the connectivity matrix
    stand_conn = standardize(sub_conn, case_masks[control])

    # Construct design matrix
    if type(control) == str:
        contrast = f'C({group}, Treatment("{control}"))'
    else:
        contrast = f'C({group}, Treatment({control}))'
        
    formula = ' + '.join((regressors, contrast))
    dmat = pat.dmatrix(formula, sub_pheno, return_type='dataframe')

    betas, pvals = glm(sub_conn, dmat, group)
    stand_betas, _ = glm(stand_conn, dmat, group)
    table = pd.DataFrame(data={'betas': betas, 'stand_betas': stand_betas, 'pvals': pvals})

    return table