import pandas as pd
import bokeh
import bebi103
from irms.constants import CONTROL


def load_data_frame_from_excel_file(input_file_path, excel_tab="d2H_data"):
    original_df = pd.read_excel(input_file_path, excel_tab)
    return original_df


def clean_up_data(original_df):
    df = original_df.rename(
        columns={
            "Identifier 1": "sample_id",
            "Rt": "relative_time",
            "Area 2": "area_2",
            "Area 3": "area_3",
            "Ampl  2": "amplitude_2",
            "Ampl  3": "amplitude_3",
            "d 3H2/2H2": "d_3H2/2H2",
            "R 3H2/2H2": "3H2_2H2_ratio",
            "d 2H/1H": "d_2H/1H",
        }
    )
    df.drop(["Rt.1"], axis=1, inplace=True)
    df["bio_replicate"] = df["sample_id"].apply(create_bio_replicate_col)
    df["TENTATIVE_COMPOUND"] = None
    df = df.groupby("Analysis", sort="relative_tine", group_keys=False).apply(
        create_original_peak_cols
    )
    df["bio_replicate_id"] = df.apply(create_short_name_col, axis=1)
    df.drop(df.index[df["bio_replicate"] == CONTROL], inplace=True)

    df = update_color_column(df, "Analysis")
    df["index_col"] = range(0, df.shape[0])

    return df


# helper functions


def create_bio_replicate_col(sample_id):
    if sample_id == "F8":
        return CONTROL
    else:
        return sample_id[:-1]


def create_original_peak_cols(grp):
    grp["original_peak_order"] = range(0, len(grp))
    grp["original_peak_count"] = len(grp)
    return grp


def create_short_name_col(row):
    if row["bio_replicate"] == CONTROL:
        return 100
    else:
        return int(row["bio_replicate"][-2:])


def update_color_column(df, column_name="Analysis"):
    df["color"] = bebi103.viz.q_to_color(
        df[column_name], colors=bokeh.palettes.Category20_10
    )
    return df
