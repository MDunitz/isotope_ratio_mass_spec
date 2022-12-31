import numpy as np
from irms.constants import amino_acids, CONTROL, NAPOI


def label_H2_and_alanine(df):
    df = handle_first_2_h2(df)
    df = handle_last_2_h2(df)
    df = handle_alanine(df)
    df = df.groupby("Analysis", group_keys=False).apply(
        create_time_relative_to_alanine_col
    )
    return df


def label_starting_H2s(analysis_grp):
    grp = analysis_grp[
        (analysis_grp["bio_replicate_id"] != CONTROL)
        & (analysis_grp["original_peak_order"] <= 4)
    ]
    h2_count = 0
    potential_h2_peak_count = len(grp[grp["amplitude_2"] > 3000])
    if potential_h2_peak_count != 2:
        analysis_id = grp["Analysis"].values[0]
        amplitudes = grp["amplitude_2"].values
        message = f"Anlaysis Group {analysis_id} has weird looking data for amplitude_2 for the first five peaks {amplitudes}"  # NOQA E501
        raise Exception(message)
    for i in range(5):
        row = grp.iloc[i]
        row_id = row["index_col"]
        if h2_count == 2:
            return analysis_grp
        if row["amplitude_2"] > 3000:
            analysis_grp.at[row_id, "TENTATIVE_COMPOUND"] = "H2"
            h2_count += 1
        else:
            analysis_grp.at[row_id, "TENTATIVE_COMPOUND"] = NAPOI


def label_final_H2s(analysis_grp):
    grp = analysis_grp[
        (analysis_grp["bio_replicate_id"] != CONTROL)
        & (
            analysis_grp["original_peak_order"]
            >= analysis_grp["original_peak_count"] - 5
        )
    ]
    h2_count = 0
    potential_h2_peak_count = len(grp[grp["amplitude_2"] > 3000])
    if potential_h2_peak_count != 2:
        analysis_id = grp["Analysis"].values[0]
        amplitudes = grp["amplitude_2"].values
        message = f"Anlaysis Group {analysis_id} has weird looking data for amplitude_2 for the first five peaks {amplitudes}"
        raise Exception(message)
    for i in range(4, 0, -1):
        row = grp.iloc[i]
        row_id = row["index_col"]
        if h2_count == 2:
            return analysis_grp
        if row["amplitude_2"] > 3000:
            analysis_grp.at[row_id, "TENTATIVE_COMPOUND"] = "H2"
            h2_count += 1
        else:
            analysis_grp.at[row_id, "TENTATIVE_COMPOUND"] = NAPOI


def label_alanine(grp):
    two_peaks = grp[(grp["peak_order"] == 2) | (grp["peak_order"] == 3)]
    potential_alanine_peak_count = len(two_peaks[two_peaks["amplitude_2"] > 900])
    if potential_alanine_peak_count < 1:
        analysis_id = two_peaks["Analysis"].values[0]
        amplitudes = two_peaks["amplitude_2"].values
        message = f"Anlaysis Group {analysis_id} has weird looking data for amplitude_2 for the first two peaks after H2 {amplitudes}"
        raise Exception(message)
    for i in range(2):
        row = two_peaks.iloc[i]
        row_id = row["index_col"]
        if row["amplitude_2"] > 900:
            grp.at[row_id, "TENTATIVE_COMPOUND"] = amino_acids["ALA"]["full_name"]
            return grp
        else:
            grp.at[row_id, "TENTATIVE_COMPOUND"] = NAPOI


def handle_alanine(df):
    df.set_index(keys="index_col", drop=False, inplace=True)
    df = df.groupby("Analysis", sort="peak_order", group_keys=False).apply(
        label_alanine
    )
    df = update_peak_order(df)
    df.reset_index(drop=True, inplace=True)
    return df


def handle_first_2_h2(df):
    df.set_index(keys="index_col", drop=False, inplace=True)
    df = df.groupby("Analysis", sort="original_peak_order", group_keys=False).apply(
        label_starting_H2s
    )
    df = update_peak_order(df)
    df.reset_index(drop=True, inplace=True)
    return df


def handle_last_2_h2(df):
    df.set_index(keys="index_col", drop=False, inplace=True)
    df = df.groupby("Analysis", sort="original_peak_order", group_keys=False).apply(
        label_final_H2s
    )
    df = update_peak_order(df)
    df.reset_index(drop=True, inplace=True)
    return df


def set_ammino_acid(df, amino_acid, time_variance, min_amplitude=0):
    source_of_truth = df[df["Analysis"] == 4167][
        [
            "TIME_RELATIVE_TO_ALANINE",
            "compound",
            "amplitude_2",
        ]
    ]
    aa_time = int(
        source_of_truth[
            source_of_truth["compound"] == amino_acids[amino_acid]["short_name"]
        ]["TIME_RELATIVE_TO_ALANINE"].values[0]
    )
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


def set_methionine(df, time_variance=7, verbose=False):
    df = df.groupby("Analysis", group_keys=False).apply(
        create_time_relative_to_glutamic_acid
    )
    df = df.groupby("Analysis", group_keys=False).apply(
        create_time_relative_to_phenylalanine
    )
    source_of_truth = df[df["Analysis"] == 4167][
        [
            "TIME_RELATIVE_TO_ALANINE",
            "TIME_RELATIVE_TO_GLUTAMIC_ACID",
            "TIME_RELATIVE_TO_PHENYLALANINE",
            "compound",
            "amplitude_2",
            "TENTATIVE_COMPOUND",
        ]
    ]
    methionine_time = int(
        source_of_truth[source_of_truth["compound"] == "Met"][
            "TIME_RELATIVE_TO_ALANINE"
        ].values[0]
    )
    methionine_time_to_GLU = int(
        source_of_truth[source_of_truth["compound"] == "Met"][
            "TIME_RELATIVE_TO_GLUTAMIC_ACID"
        ].values[0]
    )
    methionine_time_to_PHE = int(
        source_of_truth[source_of_truth["compound"] == "Met"][
            "TIME_RELATIVE_TO_PHENYLALANINE"
        ].values[0]
    )
    possible_methionine_by_alanine = df[
        (df["TIME_RELATIVE_TO_ALANINE"] <= methionine_time + time_variance)
        & (df["TIME_RELATIVE_TO_ALANINE"] >= methionine_time - time_variance)
        & (df["TENTATIVE_COMPOUND"].isnull())
    ]
    possible_methionine_by_glu = df[
        (df["TIME_RELATIVE_TO_GLUTAMIC_ACID"] <= methionine_time_to_GLU + time_variance)
        & (
            df["TIME_RELATIVE_TO_GLUTAMIC_ACID"]
            >= methionine_time_to_GLU - time_variance
        )
        & (df["TENTATIVE_COMPOUND"].isnull())
    ]
    possible_methionine_by_phe = df[
        (df["TIME_RELATIVE_TO_PHENYLALANINE"] <= methionine_time_to_PHE + time_variance)
        & (
            df["TIME_RELATIVE_TO_PHENYLALANINE"]
            >= methionine_time_to_PHE - time_variance
        )
        & (df["TENTATIVE_COMPOUND"].isnull())
    ]

    df.loc[
        possible_methionine_by_alanine.index, ["TENTATIVE_COMPOUND_BASED_ON_ALANINE"]
    ] = amino_acids["MET"]["full_name"]
    df.loc[
        possible_methionine_by_glu.index, ["TENTATIVE_COMPOUND_BASED_ON_GLU"]
    ] = amino_acids["MET"]["full_name"]
    df.loc[
        possible_methionine_by_phe.index, ["TENTATIVE_COMPOUND_BASED_ON_PHE"]
    ] = amino_acids["MET"]["full_name"]

    ala_set = set(possible_methionine_by_alanine.index)
    phe_set = set(possible_methionine_by_phe.index)
    glu_set = set(possible_methionine_by_glu.index)
    potential_indexes = np.array(
        list(ala_set.intersection(phe_set).intersection(glu_set))
    )
    possible_mets = df.loc[potential_indexes]
    if verbose:
        print(f"total sample count: {len(df['sample_id'].unique())}")
        print(f"time variance: {time_variance}")
        print(
            f"met by ala: {len(possible_methionine_by_alanine['sample_id'].values)}",
            f"unique sample ids: {len(possible_methionine_by_alanine['sample_id'].unique())}",
        )
        print(
            f"met by glu: {len(possible_methionine_by_glu['sample_id'].values)}",
            f"unique sample ids: {len(possible_methionine_by_glu['sample_id'].unique())}",
        )
        print(
            f"met by phe: {len(possible_methionine_by_phe['sample_id'].values)}",
            f"unique sample ids: {len(possible_methionine_by_phe['sample_id'].unique())}",
        )

        print(
            f"index intersection len: {len(ala_set.intersection(phe_set).intersection(glu_set))}"
        )
        print(f"index union len: {len(ala_set.union(phe_set).union(glu_set))}")
        print(
            f"ala mean: \n{possible_methionine_by_alanine[['amplitude_2', 'area_2']].mean()}"
        )
        print(
            f"glu mean: \n{possible_methionine_by_glu[['amplitude_2', 'area_2']].mean()}"
        )
        print(
            f"phe mean: \n{possible_methionine_by_phe[['amplitude_2', 'area_2']].mean()}"
        )
        print(
            f"unique samples in met intersection: {len(possible_mets['sample_id'].unique())}"
        )
        print(possible_mets.shape[0])
    assert len(possible_mets["sample_id"].unique()) == len(possible_mets["sample_id"])
    print(f"There are {len(possible_mets)} possible Mets")
    df.loc[possible_mets.index, ["TENTATIVE_COMPOUND"]] = amino_acids["MET"][
        "full_name"
    ]
    return df


def label_required_aas(df):
    df = set_ammino_acid(df, "VAL", 6)
    df = set_ammino_acid(df, "LEU", 3)
    df = set_ammino_acid(df, "ILE", 5, 600)
    df = set_ammino_acid(df, "PRO", 5)
    df = set_ammino_acid(df, "ASP", 9, 700)
    df = set_ammino_acid(df, "PHE", 6, 500)
    return df


def label_non_required_aas(df):
    df = set_ammino_acid(df, "GLU", 6)
    df = set_ammino_acid(df, "LYS", 4)
    df = set_ammino_acid(df, "SER", 8)
    df = set_ammino_acid(df, "TYR", 8)
    df = set_methionine(df, 7)
    return df


# helper functions


def create_time_relative_to_alanine_col(grp):
    alanine_time = int(
        grp[grp["TENTATIVE_COMPOUND"] == amino_acids["ALA"]["full_name"]][
            "relative_time"
        ]
    )
    grp["TIME_RELATIVE_TO_ALANINE"] = grp["relative_time"] - alanine_time
    return grp


def create_time_relative_to_glutamic_acid(grp):
    glutamic_acid_time = int(
        grp[grp["TENTATIVE_COMPOUND"] == amino_acids["GLU"]["full_name"]][
            "relative_time"
        ]
    )
    grp["TIME_RELATIVE_TO_GLUTAMIC_ACID"] = grp["relative_time"] - glutamic_acid_time
    return grp


def create_time_relative_to_phenylalanine(grp):
    phenylalanine_time = int(
        grp[grp["TENTATIVE_COMPOUND"] == amino_acids["PHE"]["full_name"]][
            "relative_time"
        ]
    )
    grp["TIME_RELATIVE_TO_PHENYLALANINE"] = grp["relative_time"] - phenylalanine_time
    return grp


def reset_peak_order_col(grp):
    grp["peak_order"] = range(0, len(grp))
    grp["peak_count"] = len(grp)
    return grp


def update_peak_order(df):
    print("Updating Peak Order", df.shape)
    df.drop(df.index[df["TENTATIVE_COMPOUND"] == NAPOI], inplace=True)
    print("after", df.shape)
    df = df.groupby("Analysis", sort="relative_tine", group_keys=False).apply(
        reset_peak_order_col
    )
    return df
