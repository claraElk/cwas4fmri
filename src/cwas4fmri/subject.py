import os
import json
import warnings
from tqdm import tqdm
import numpy as np

def find_valid_subjects(bids_dir, pheno, session, connectome_t, run, task, atlas, feature, out_p):    
    print("⏳ Identify subjects connectivity matrix ...")
    print("This will take a moment, please do not interupt the process ...\n")
    
    # Find all subjects in phenotype file
    all_subjects = set(pheno['participant_id'])
    
    # Find subjects processed by HALFpipe
    processed_subjects = set()

    for _, row in tqdm(pheno.iterrows()):
        file_path = os.path.join(bids_dir, connectome_t.format(row['participant_id'], session, 
                                                               row['participant_id'], session, 
                                                               task, run, atlas, feature)
                                                                )
        if os.path.exists(file_path):
            processed_subjects.add(row['participant_id'])
        
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
    report_file = os.path.join(out_p, 'cwas_report.json')
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=4)
    
    # Print results
    print("\n=== Summary subjects processing ===")
    print(f"Total subjects: {len(all_subjects)}")
    print(f"Processed subjects: {len(processed_subjects)}")
    print(f"Unprocessed subjects: {len(unprocessed_subjects)}")
    
    if unprocessed_subjects:
        print("\n❗️ Subjects not processed:")
        for subject in sorted(unprocessed_subjects):
            print(f"- {subject}")
    else:
        print("✅ All subjects have been processed!")

    print(f"\n✅ Information saved in:", report_file)
    
    return pheno[pheno['participant_id'].isin(processed_subjects)].reset_index(drop=True)


def find_subset(pheno, column, cases=None):
    subset_mask = np.array(~pheno[column].isnull())
    
    if cases is not None and not not cases:
        all_cases = pheno.loc[subset_mask][column].unique()

        try:
            case_available = np.array([True if case in all_cases else False for case in cases])

        except TypeError as e:
            raise Exception(f'❌ the attribute "cases" needs to be iterable but is: {type(cases)}')
        
        if not all(case_available):
            if not any(case_available):
                raise Exception(f'❌ none of the requested cases {cases} of "{column}" are available')
            else:
                warnings.warn(
                    f'\n❌ not all requested cases of "{column}" are available: {list(zip(cases, case_available))}',
                    RuntimeWarning)
                
        case_masks = np.array([pheno[column] == case for case in cases])
        subset_mask = np.any(case_masks, 0)
            
        # Return the masked instances of the requested cases
        cases_dict = {case: case_masks[idx][subset_mask] for idx, case in enumerate(cases)}
        return subset_mask, cases_dict
    else:
        return subset_mask