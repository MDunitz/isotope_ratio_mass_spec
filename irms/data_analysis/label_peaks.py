from irms.constants import amino_acids, name_to_short_name
from irms.data_analysis import first_pass



def get_aa_mean_relative_time(df):
    aa_peak_summary_data = {}
    for group, data in df.groupby("TENTATIVE_COMPOUND"):
        short_name = name_to_short_name[group]
        aa_peak_summary_data[short_name] = {}
        aa_peak_summary_data[short_name]['mean_alanine_relative_time'] = data["TIME_RELATIVE_TO_ALANINE"].mean().round()
    aa_peak_summary_data.pop("H2")
    aa_peak_summary_data.pop("ALA")
    return aa_peak_summary_data




def get_time_variances(aa_peak_summaries, df):
    for amino_acid, aa_data in aa_peak_summaries.items():
        time_variance = find_time_variance(df, amino_acid, aa_data['mean_alanine_relative_time'])
        aa_peak_summaries[amino_acid]['time_variance'] = time_variance
    return aa_peak_summaries


def set_ammino_acid(df, aa_time, amino_acid, time_variance, min_amplitude=0):
    possible_aas = df[
        (df["TIME_RELATIVE_TO_ALANINE"] <= aa_time + time_variance)
        & (df["TIME_RELATIVE_TO_ALANINE"] >= aa_time - time_variance)
        & (df["TENTATIVE_COMPOUND"].isnull())
        & (df["amplitude_2"] >= min_amplitude)
    ]
    print(f"There are {len(possible_aas)} possible {amino_acid}s")
    for group, data in possible_aas.groupby("Analysis"):
        if len(data) == 1:
            pass
        else:
            raise Exception(f"More than one {amino_acid} found in {group}")
    df.loc[possible_aas.index, ["TENTATIVE_COMPOUND"]] = amino_acids[amino_acid][
        "full_name"
    ]
    return df




def label_aas(df, aa_peak_summaries):
    print(aa_peak_summaries)
    for amino_acid, data in aa_peak_summaries.items():
        min_amplitude = amino_acids[amino_acid]["min_amp"]
        print(amino_acid, data)
        df = set_ammino_acid(df, data['mean_alanine_relative_time'], amino_acid, data['time_variance'], min_amplitude)
    return df



def find_time_variance(df, amino_acid_name, alanine_relative_time):
    analysis_runs = len(df["Analysis"].unique())
    min_amplitude = amino_acids[amino_acid_name]["min_amp"]
    for i in range(10):
        analysis_ids = get_possible_rows_for_aa(df, alanine_relative_time, i, min_amplitude)['Analysis']
        if analysis_ids.is_unique:
            if len(analysis_ids) == analysis_runs:
                print(f"{i} for: {amino_acid_name}, peak_count: {len(analysis_ids)}, unique pc: {len(analysis_ids.unique())}")
                return i
        else:
            duplicates = [x[0] for x in analysis_ids.value_counts().items() if x[1] > 1]
            print(f"{i} for: {amino_acid_name}, peak_count: {len(analysis_ids)}, unique pc: {len(analysis_ids.unique())}")
            print(f"More than one {amino_acid_name} found at time variance {i} in analysis groups: {duplicates}")
            return i-1
    return i

def get_possible_rows_for_aa(df, alanine_relative_time, time_variance, min_amplitude=0):
    possible_aas = df[
        (df["TIME_RELATIVE_TO_ALANINE"] <= alanine_relative_time + time_variance)
        & (df["TIME_RELATIVE_TO_ALANINE"] >= alanine_relative_time - time_variance)
        & (df["TENTATIVE_COMPOUND"].isnull())
        & (df["amplitude_2"] >= min_amplitude)
    ]
    return possible_aas
