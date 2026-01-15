import sys
import time
import warnings
import scipy as sp
import numpy as np
import patsy as pat
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from sklearn import linear_model as sln
from sklearn import preprocessing as skp
from sklearn import model_selection as skm

def conn2mat(flat_conn, mask_2d):
    """
    Convert flattened upper-triangle vector to full square matrix.
    """
    n_rois = mask_2d.shape[0]
    mat = np.zeros((n_rois, n_rois), dtype=np.float32)
    mat[mask_2d] = flat_conn
    # Mirror to lower triangle
    mat = mat + mat.T
    return mat


def find_subset(pheno, column, cases=None):
    subset_mask = np.array(~pheno[column].isnull())
    
    if cases is not None and not not cases:
        all_cases = pheno.loc[subset_mask][column].unique()
        
        try:
            case_available = np.array([True if case in all_cases else False for case in cases])
        except TypeError as e:
            raise Exception(f'the attribute "cases" needs to be iterable but is: {type(cases)}') from e
            
        if not all(case_available):
            if not any(case_available):
                raise Exception(f'none of the requested cases of "{column}" are available')
            else:
                warnings.warn(
                    f'\nnot all requested cases of "{column}" are available: {list(zip(cases, case_available))}',
                    RuntimeWarning)
                
        case_masks = np.array([pheno[column] == case for case in cases])
        subset_mask = np.any(case_masks, 0)
        # Return the masked instances of the requested cases
        cases_dict = {case: case_masks[idx][subset_mask] for idx, case in enumerate(cases)}
        return subset_mask, cases_dict
    else:
        return subset_mask


def standardize(data, mask):
    """
    Standardize based on HC
    """
    scaler = skp.StandardScaler(with_mean=False, with_std=True)
    scaler.fit(data[mask, :])
    standardized_data = scaler.transform(data)
    return standardized_data


def find_contrast(design_matrix, contrast):
    """
    Find the contrast column
    """
    contrast_columns = [(col_id, col) for col_id, col in enumerate(design_matrix.columns) if f'{contrast}' in col]    
    if not len(contrast_columns) == 1:
        raise Exception(f'There is no single factor that matches {contrast}: {(list(design_matrix.columns))}')
    return contrast_columns


def fast_glm(data, design_matrix, contrast):
    # Does not compute p-values but operates in parallel
    contrast_id, contrast_name = find_contrast(design_matrix, contrast)[0]
    glm = sln.LinearRegression(fit_intercept=False, normalize=False, n_jobs=-2)
    res = glm.fit(design_matrix, data)
    betas = res.coef_[:, contrast_id]
    return betas


def glm(data, design_matrix, contrast):
    """
    Perform GLM
    """
    contrast_id, contrast_name = find_contrast(design_matrix, contrast)[0]
    n_data = data.shape[1]

    # Conduct the GLM
    betas = np.zeros(shape=n_data)
    pvals = np.zeros(shape=n_data)
    for conn_id in range(n_data):
        model = sm.OLS(data[:, conn_id], design_matrix, missing = 'drop')
        results = model.fit()
        betas[conn_id] = results.params.iloc[contrast_id]
        pvals[conn_id] = results.pvalues.iloc[contrast_id]


    return betas, pvals


def glm_wrap_cc(conn, pheno, group, case, control, regressors='', report=False, fast=False):
    # Make sure pheno and conn have the same number of cases
    if not conn.shape[0] == pheno.shape[0]:
        print(f'Conn ({conn.shape[0]}) and pheno ({pheno.shape[0]}) must be same number of cases')

    # Define the subset of the sample
    sub_mask, case_masks = find_subset(pheno, group, [case, control])
    sub_conn = conn[sub_mask, :]
    sub_pheno = pheno.loc[sub_mask]
    n_sub = np.sum(sub_mask)
    n_case = np.sum(case_masks[case])
    n_control = np.sum(case_masks[control])
    n_data = sub_conn.shape[1]
    
    if report:
        print(f'Selected sample based on group variable {group}.\n'
              f'cases: {case} (n={n_case})\n'
              f'controls: {control} (n={n_control})\n'
              f'original sample: n={pheno.shape[0]}; new sample: n={n_sub}\n'
              f'{n_data} data points available\n'
              f'standardized estimators are based on {group}=={control}')

    stand_conn = standardize(sub_conn, case_masks[control])
    
    # Construct design matrix
    if type(control) == str:
        contrast = f'C({group}, Treatment("{control}"))'
    else:
        contrast = f'C({group}, Treatment({control}))'
        
    formula = ' + '.join((regressors, contrast))
    dmat = pat.dmatrix(formula, sub_pheno, return_type='dataframe')

    if fast: # We don't use that for now
        betas = fast_glm(sub_conn, dmat, group)
        table = pd.DataFrame(data={'betas': betas})
    else:
        betas, pvals = glm(sub_conn, dmat, group)
        stand_betas, _ = glm(stand_conn, dmat, group)
        table = pd.DataFrame(data={'betas': betas, 'stand_betas': stand_betas, 'pvals': pvals})

    return table