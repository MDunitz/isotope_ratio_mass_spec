import numpy as np
from irms.constants import amino_acids, CONTROL, NAPOI




def set_ammino_acid(df, control, amino_acid, time_variance, min_amplitude=0):
    aa_time = int(
        control[
            control["compound"] == amino_acids[amino_acid]["short_name"]
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
    control = df[df["Analysis"] == 4167][
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
        control[control["compound"] == "Met"]["TIME_RELATIVE_TO_ALANINE"].values[0]
    )
    methionine_time_to_GLU = int(
        control[control["compound"] == "Met"]["TIME_RELATIVE_TO_GLUTAMIC_ACID"].values[
            0
        ]
    )
    methionine_time_to_PHE = int(
        control[control["compound"] == "Met"]["TIME_RELATIVE_TO_PHENYLALANINE"].values[
            0
        ]
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


def label_required_aas(df, control=None):
    if control is None:
        control = df[df["Analysis"] == 4167][
            [
                "TIME_RELATIVE_TO_ALANINE",
                "compound",
                "amplitude_2",
            ]
        ]
    df = set_ammino_acid(df, control, "VAL", 6)
    df = set_ammino_acid(df, control, "LEU", 3)
    df = set_ammino_acid(df, control, "ILE", 5, 600)
    df = set_ammino_acid(df, control, "PRO", 5)
    df = set_ammino_acid(df, control, "ASP", 9, 700)
    df = set_ammino_acid(df, control, "PHE", 6, 500)
    return df


def label_non_required_aas(df, control=None):
    if control is None:
        control = df[df["Analysis"] == 4167][
            [
                "TIME_RELATIVE_TO_ALANINE",
                "compound",
                "amplitude_2",
            ]
        ]
    df = set_ammino_acid(df, control, "GLU", 6)
    df = set_ammino_acid(df, control, "LYS", 4)
    df = set_ammino_acid(df, control, "SER", 8)
    df = set_ammino_acid(df, control, "TYR", 8)
    df = set_methionine(df, 7)
    return df


# helper functions


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

