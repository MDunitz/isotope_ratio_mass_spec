import time
import pandas as pd
from irms import data_wrangling
from irms.constants import amino_acids
from irms.data_analysis import label_peaks
from irms.data_analysis import first_pass
from irms.data_analysis import label_h2_and_alanine

def label_aa_peaks_based_on_manual_analysis(path_to_data, write=False, output_file=None):
    original_df = data_wrangling.load_data_frame_from_excel_file(path_to_data)
    df = data_wrangling.clean_up_data(original_df)
    df = label_h2_and_alanine.label_H2_and_alanine(df)
    df = first_pass.label_required_aas(df)
    df = first_pass.label_non_required_aas(df)
    # TODO clean up df (remove rando columns before writing to file)
    if write:
        if output_file is None:
            time = time.time().split('.')[0]
            output_file = f"you_forgot_to_add_a_filename-{time}.csv"
        df.to_csv(output_file, index=False)
    return df


def label_aa_peaks_based_on_avg_time_to_ALA(path_to_original_data, write=False, output_file=None, labeled_df=None, path_to_labeled_data=None):
    original_df = data_wrangling.load_data_frame_from_excel_file(path_to_original_data)
    clean_df = data_wrangling.clean_up_data(original_df)
    df = label_h2_and_alanine.label_H2_and_alanine(clean_df)
    # either pass in the labeled df, or a path to a file with the labeled df or regenerate it from the original data
    if labeled_df is None:
        if path_to_labeled_data:
            labeled_df = pd.read_csv(path_to_labeled_data)
        else:
            labeled_df = label_aa_peaks_based_on_manual_analysis(path_to_data=path_to_original_data)
            
    aa_mean_relative_time = label_peaks.get_aa_mean_relative_time(labeled_df)
    aa_peak_summaries = label_peaks.get_time_variances(aa_mean_relative_time, df)
    df = label_peaks.label_aas(df, aa_peak_summaries)
    if write:
        if output_file is None:
            time = time.time().split('.')[0]
            output_file = f"you_forgot_to_add_a_filename-{time}.csv"
        df.to_csv(output_file, index=False)
    return df

    