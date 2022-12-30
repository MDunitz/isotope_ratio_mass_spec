import pandas as pd
from irms import data_wrangling
from irms.data_analysis import label_peaks
from irms.constants import *


def label_aa_peaks(path_to_data, output_file):
    original_df = data_wrangling.load_data_frame_from_excel_file(path_to_data)
    df = data_wrangling.clean_up_data(original_df)
    df = label_peaks.label_H2_and_alanine(df)
    df = label_peaks.label_required_aas(df)
    df = label_peaks.label_non_required_aas(df)
    ## Todo clean up df (remove rando columns before writing to file)
    df.to_csv(output_file, index=False)
    return df
