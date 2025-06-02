# In python by default
import warnings
from pathlib import Path

# Required installation
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)


def create_output(out_p) : 
    # Create output directory if it doesn't exist
    out_p = Path(out_p)
    out_p.mkdir(parents=True, exist_ok=True)

def read_phenotype_file(path):
    ext = Path(path).suffix.lower()
    if ext not in ['.tsv', '.csv', '.xlsx']:
        raise ValueError(f"Unsupported extension: {ext}")
    try:
        if ext == '.tsv':
            return pd.read_csv(path, sep='\t')
        elif ext == '.csv':
            return pd.read_csv(path)
        else:
            return pd.read_excel(path)
    except Exception as e:
        raise ValueError(f"Error reading file: {e}")


def validate_columns(df, col_map):
    missing = [c for c in col_map.keys() if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    df = df.rename(columns=col_map)
    return df


def validate_subject_ids(df, subject_file_path):
    with open(subject_file_path) as f:
        subjects = [line.strip() for line in f if line.strip()]
    pattern = r'^sub-(?:' + '|'.join(subjects) + r')$'
    invalid = df[~df['participant_id'].str.match(pattern)]
    if not invalid.empty:
        warnings.warn(f"Found invalid subject IDs: {invalid['participant_id'].tolist()[:5]}")
    return df


def encode_diagnosis(df, case, control):
    allowed = [case, control]
    invalid = df[~df["diagnosis"].isin(allowed)]
    if not invalid.empty:
        warnings.warn(f"Unexpected diagnosis values: {invalid['diagnosis'].unique()}")
    df = df[df["diagnosis"].isin(allowed)]
    return df.replace({'diagnosis': {case: 0, control: 1}})


def encode_sex(df):
    # TODO: check if only 2 unique values are present
    sex_variable = df["sex"].unique()
    print(sex_variable)

    return df.replace({'sex': {sex_variable[0]: 0, sex_variable[1]: 1}})


def warn_on_invalid_age(df):
    invalid = df[~df["age"].apply(lambda x: isinstance(x, (int, float)) and x > 0)]
    if not invalid.empty:
        warnings.warn(f"Possible invalid age values in rows: {invalid.index.tolist()[:5]}")

def load_phenotype(
    phenotype_file_path, diagnosis_col, subject_col, age_col, sex_col,
    scanner_col=False, sequence_col=False, medication_col=False,
    case_name="case", control_name="control",
    subject_file_path="subject-list.txt"
):
    print(f"Reading phenotype file: {phenotype_file_path}")
    df = read_phenotype_file(phenotype_file_path)

    # TODO : once approved, delete renaming of columns
    col_map = {
        subject_col: "participant_id",
        diagnosis_col: "diagnosis",
        age_col: "age",
        sex_col: "sex"
    }

    # Hard coded columns
    if scanner_col: 
        col_map["scanner"] = "scanner"
    if sequence_col: 
        col_map["sequence"] = "sequence"
    if medication_col: 
        col_map["medication"] = "medication"

    df = validate_columns(df, col_map)

    if "scanner" in df.columns:
        print("Scanners found:", df["scanner"].unique())
    if "sequence" in df.columns:
        print("Sequences found:", df["sequence"].unique())
    if "medication" in df.columns:
        print("Medications found:", df["medication"].unique())

    df = validate_subject_ids(df, subject_file_path)
    df = encode_diagnosis(df, case_name, control_name)
    df = encode_sex(df)
    warn_on_invalid_age(df)

    return df