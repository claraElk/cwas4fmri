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
        raise ValueError(f"❌ Unsupported extension: {ext}. Expected .tsv, .csv, or .xlsx\n")
    try:
        if ext == '.tsv':
            return pd.read_csv(path, sep='\t')
        elif ext == '.csv':
            return pd.read_csv(path)
        else:
            return pd.read_excel(path)
    except Exception as e:
        raise ValueError(f"❌ Error reading file: {e}. Make sure the file exists and is in the correct format.\n")


def validate_columns(df, col_map):
    missing = [c for c in col_map.keys() if c not in df.columns]
    if missing:
        raise ValueError(f"❌ Missing required columns: {missing}\n")
    df = df.rename(columns=col_map)
    return df


def validate_subject_ids(df, subject_file_path):
    with open(subject_file_path) as f:
        subjects = [line.strip() for line in f if line.strip()]
    pattern = r'^sub-(?:' + '|'.join(subjects) + r')$'
    invalid = df[~df['participant_id'].str.match(pattern)]
    if not invalid.empty:
        warnings.warn(f"❗️ Found invalid subject IDs: {invalid['participant_id'].tolist()[:5]} \n")
    return df


def encode_diagnosis(df, case, control):
    allowed = [case, control]
    print(f"Encoding diagnosis: {case} as 0, {control} as 1")
    invalid = df[~df["diagnosis"].isin(allowed)]
    if not invalid.empty:
        warnings.warn(f"❗️ Unexpected diagnosis values: {invalid['diagnosis'].unique()}\n")
    df = df[df["diagnosis"].isin(allowed)]
    return df.replace({'diagnosis': {case: 0, control: 1}})


def encode_sex(df):
    sex_variable = df["sex"].unique()
    return df.replace({'sex': {sex_variable[0]: 0, sex_variable[1]: 1}})


def warn_on_invalid_age(df):
    invalid = df[~df["age"].apply(lambda x: isinstance(x, (int, float)) and x > 0)]
    if not invalid.empty:
        warnings.warn(f"❗️ Possible invalid age values in rows: {invalid.index.tolist()[:5]}\n")

def load_phenotype(
    phenotype_file_path, diagnosis_col, subject_col, age_col, sex_col,
    scanner_col=False, sequence=False, medication=False,
    case_name="case", control_name="control"
):
    print(f"⏳ Reading phenotype file: {phenotype_file_path}")
    df = read_phenotype_file(phenotype_file_path)

    col_map = {
        subject_col: "participant_id",
        diagnosis_col: "diagnosis",
        age_col: "age",
        sex_col: "sex"
    }

    # Hard coded columns
    if scanner_col: 
        col_map["scanner"] = "scanner"
    if sequence: 
        col_map["sequence"] = "sequence"
    if medication: 
        col_map["medication"] = "medication"

    df = validate_columns(df, col_map)

    if scanner_col:
        print("Scanners found:", df["scanner"].unique())
    if sequence:
        print("Sequences found:", df["sequence"].unique())
    if medication:
        print("Medications found:", df["medication"].unique())

    df = encode_diagnosis(df, case_name, control_name)
    df = encode_sex(df)
    warn_on_invalid_age(df)

    return df