import glob
import os
import json
import warnings
from tqdm import tqdm
import numpy as np

def find_valid_subjects(pheno, derivatives_p, confounds_p, out_p):    
    print("Identify subjects connectivity matrix ...")
    print("This will take a moment, please do not interupt the process ...")
    
    # Find all subjects in phenotype file
    all_subjects = set(pheno['participant_id'])
    
    # Find subjects processed by HALFpipe
    processed_subjects = set()
    for rid, row in tqdm(pheno.iterrows()):
        file_path = glob.glob(os.path.join(derivatives_p, confounds_p.format(row['participant_id'], row['participant_id'])),
                             recursive=True)
        
        if any(os.path.exists(path) for path in file_path):
            processed_subjects.add(row['participant_id'])
    
    print("print file_path", file_path)
    
    # Find unprocessed subjects
    unprocessed_subjects = all_subjects - processed_subjects
    
    # Create report
    report = {
        "halfpipe_unprocessed": {
            "total_subjects": len(all_subjects),
            "processed_subjects": len(processed_subjects),
            "unprocessed_subjects": len(unprocessed_subjects),
            "unprocessed_subject_list": list(unprocessed_subjects)
        }
    }
    
    # Create results directory if it doesn't exist
    os.makedirs(out_p, exist_ok=True)
    
    # Save report to JSON
    report_file = os.path.join(out_p, 'cwas_subject_report.json')
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=4)
    
    # Print results
    print("\n=== Summary HALFpipe processing ===")
    print(f"Total subjects: {len(all_subjects)}")
    print(f"Processed subjects: {len(processed_subjects)}")
    print(f"Unprocessed subjects: {len(unprocessed_subjects)}")
    
    if unprocessed_subjects:
        print("\nSubjects not processed by HALFpipe:")
        for subject in sorted(unprocessed_subjects):
            print(f"- {subject}")
    else:
        print("All subjects have been processed by HALFpipe!")

    print(f"Information saved in:", report_file)
    
    return pheno[pheno['participant_id'].isin(processed_subjects)].reset_index(drop=True)


def find_subset(pheno, column, cases=None):
    subset_mask = np.array(~pheno[column].isnull())
    
    if cases is not None and not not cases:
        print("in cases")
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
