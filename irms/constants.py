CONTROL = "CONTROL"
REQUIRED = "REQUIRED"
NAPOI = "NOT A PEAK OF INTEREST"

amino_acids = {
    "ALA": {"full_name": "Alanine", "short_name": "Ala", "likelihood": 1, "min_amp": 0},
    "ARG": {"full_name": "Arginine", "short_name": "Arg", "likelihood": 0, "min_amp": 0},
    "ASN": {"full_name": "Asparagine", "short_name": "Asn", "likelihood": 0, "min_amp": 0},
    "ASP": {"full_name": "Aspartic acid", "short_name": "Asp", "likelihood": 1, "min_amp": 700},
    "CYS": {"full_name": "Cysteine", "short_name": "Cys", "likelihood": 0, "min_amp": 0},
    "GLN": {"full_name": "Glutamine", "short_name": "Gln", "likelihood": 0, "min_amp": 0},
    "GLU": {"full_name": "Glutamic acid", "short_name": "Glu", "likelihood": 0.95, "min_amp": 0},
    "GLY": {"full_name": "Glycine", "short_name": "Gly", "likelihood": 0, "min_amp": 0},
    "HIS": {"full_name": "Histidine", "short_name": "His", "likelihood": 0, "min_amp": 0},
    "ILE": {"full_name": "Isoleucine", "short_name": "Ile", "likelihood": 1, "min_amp": 600},
    "LEU": {"full_name": "Leucine", "short_name": "Leu", "likelihood": 1, "min_amp": 0},
    "LYS": {"full_name": "Lysine", "short_name": "Lys", "likelihood": 0.92, "min_amp": 0},
    "MET": {"full_name": "Methionine", "short_name": "Met", "likelihood": 0.95, "min_amp": 0},
    "PHE": {"full_name": "Phenylalanine", "short_name": "Phe", "likelihood": 1, "min_amp": 500},
    "PRO": {"full_name": "Proline", "short_name": "Pro", "likelihood": 1, "min_amp": 0},
    "SER": {"full_name": "Serine", "short_name": "Ser", "likelihood": 0.93, "min_amp": 0},
    "THR": {"full_name": "Threonine", "short_name": "Thr", "likelihood": 0, "min_amp": 0},
    "TRP": {"full_name": "Tryptophan", "short_name": "Trp", "likelihood": 0, "min_amp": 0},
    "TYR": {"full_name": "Tyrosine", "short_name": "Tyr", "likelihood": 0.93, "min_amp": 0},
    "VAL": {"full_name": "Valine", "short_name": "Val", "likelihood": 1, "min_amp": 0},
}

mean_time_relative_to_alanine = {
    "Aspartic acid": 674.0,
    "Glutamic acid": 898.0,
    "Isoleucine": 461.0,
    "Leucine": 435.0,
    "Lysine": 1597.0,
    "Methionine": 926.0,
    "Phenylalanine": 1136.0,
    "Proline": 534.0,
    "Serine": 792.0,
    "Tyrosine": 1831.0,
    "Valine": 270.0,
}

name_to_short_name = {
    "Alanine": "ALA",
    "Aspartic acid": "ASP",
    "Glutamic acid": "GLU",
    "Isoleucine": "ILE",
    "Leucine": "LEU",
    "Lysine": "LYS",
    "Methionine": "MET",
    "Phenylalanine": "PHE",
    "Proline": "PRO",
    "Serine": "SER",
    "Tyrosine": "TYR",
    "Valine": "VAL",
    "H2": "H2"
}