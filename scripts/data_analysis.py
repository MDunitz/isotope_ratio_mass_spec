from constants import *
import pandas as pd
peaks_liklihood = {H2: 4, ALA: 1, VAL: 1, LEU: 1, ILE: 1, PRO:1, ASP:1, PHE:1,  GLU: 0.9,  MET:0.9,  SER: .75, LYS: 0.67, TYR: 0.67}


def tentatively_label_H2_and_alanine(df):
    df = handle_first_2_h2(df)
    df = handle_last_2_h2(df)
    df = handle_alanine(df)
    df = df.groupby('Analysis', group_keys=False).apply(create_time_relative_to_alanine_col)
    return df


def create_time_relative_to_alanine_col(grp):
    alanine_time = int(grp[grp['TENTATIVE_COMPOUND'] == ALA]['relative_time'])
    grp['TIME_RELATIVE_TO_ALANINE'] = grp['relative_time'] -alanine_time
    return grp

def create_time_relative_to_glutamic_acid(grp):
    glutamic_acid_time = int(grp[grp['TENTATIVE_COMPOUND'] == GLU]['relative_time'])
    grp['TIME_RELATIVE_TO_GLUTAMIC_ACID'] = grp['relative_time'] -glutamic_acid_time
    return grp 

def create_time_relative_to_phenylalanine(grp):
    phenylalanine_time = int(grp[grp['TENTATIVE_COMPOUND'] == PHE]['relative_time'])
    grp['TIME_RELATIVE_TO_PHENYLALANINE'] = grp['relative_time'] -phenylalanine_time
    return grp 

def reset_peak_order_col(grp):
    grp['peak_order'] = range(0, len(grp))
    grp['peak_count'] = len(grp)
    return grp

def update_peak_order(df):
    print("Updating Peak Order", df.shape)
    df.drop(df.index[df['TENTATIVE_COMPOUND'] == NAPOI], inplace=True)
    print("after", df.shape)
    df = df.groupby('Analysis', sort='relative_tine', group_keys=False).apply(reset_peak_order_col)
    return df


def label_starting_H2s(analysis_grp):
    grp = analysis_grp[(analysis_grp['bio_replicate_id'] != CONTROL) & (analysis_grp['original_peak_order'] <= 4)] 
    h2_count = 0
    potential_h2_peak_count = len(grp[grp['amplitude_2'] > 3000])
    if potential_h2_peak_count !=2:
        analysis_id = grp['Analysis'].values[0]
        amplitudes = grp['amplitude_2'].values
        raise Exception(
            f"Anlaysis Group {analysis_id} has weird looking data for amplitude_2 for the first five peaks {amplitudes}")
    for i in range(5):
        row = grp.iloc[i]
        row_id = row['index_col']
        if h2_count == 2:
            return analysis_grp
        if row['amplitude_2'] >3000:
            analysis_grp.at[row_id, "TENTATIVE_COMPOUND"] = "H2"
            h2_count +=1
        else:
            analysis_grp.at[row_id, "TENTATIVE_COMPOUND"] = NAPOI
        

def label_final_H2s(analysis_grp):
    grp = analysis_grp[(analysis_grp['bio_replicate_id'] != CONTROL) & (analysis_grp['original_peak_order'] >= analysis_grp['original_peak_count'] - 5)] 
    h2_count = 0
    potential_h2_peak_count = len(grp[grp['amplitude_2'] > 3000])
    if potential_h2_peak_count !=2:
        analysis_id = grp['Analysis'].values[0]
        amplitudes = grp['amplitude_2'].values
        raise Exception(
            f"Anlaysis Group {analysis_id} has weird looking data for amplitude_2 for the first five peaks {amplitudes}")
    for i in range(4, 0, -1):
        row = grp.iloc[i]
        row_id = row['index_col']
        if h2_count == 2:
            return analysis_grp
        if row['amplitude_2'] >3000:
            analysis_grp.at[row_id, "TENTATIVE_COMPOUND"] = "H2"
            h2_count +=1
        else:
            analysis_grp.at[row_id, "TENTATIVE_COMPOUND"] = NAPOI

def label_alanine(grp):
    two_peaks = grp[(grp['peak_order']==2) | (grp['peak_order']==3)]
    potential_alanine_peak_count = len(two_peaks[two_peaks['amplitude_2'] > 900])
    if potential_alanine_peak_count < 1:
        analysis_id = two_peaks['Analysis'].values[0]
        amplitudes = two_peaks['amplitude_2'].values
        raise Exception(
            f"Anlaysis Group {analysis_id} has weird looking data for amplitude_2 for the first two peaks after H2 {amplitudes}")
    for i in range(2):
        row = two_peaks.iloc[i]
        row_id = row['index_col']
        if row['amplitude_2'] > 900:
            grp.at[row_id, "TENTATIVE_COMPOUND"] = ALA
            return grp
        else:
            grp.at[row_id, "TENTATIVE_COMPOUND"] = NAPOI
        
def handle_alanine(df):
    df.set_index(keys='index_col', drop=False, inplace=True)
    df = df.groupby('Analysis', sort='peak_order', group_keys=False).apply(label_alanine)
    df = update_peak_order(df)
    df.reset_index(drop=True, inplace=True)
    return df

def handle_first_2_h2(df):
    df.set_index(keys='index_col', drop=False, inplace=True)
    df = df.groupby('Analysis', sort='original_peak_order', group_keys=False).apply(label_starting_H2s)
    df = update_peak_order(df)
    df.reset_index(drop=True, inplace=True)
    return df

def handle_last_2_h2(df):
    df.set_index(keys='index_col', drop=False, inplace=True)
    df = df.groupby('Analysis', sort='original_peak_order', group_keys=False).apply(label_final_H2s)
    df = update_peak_order(df)
    df.reset_index(drop=True, inplace=True)
    return df


## TODO move to sep file


def set_valines(df):
    TIME_VARIANCE = 6
    source_of_truth = df[df['Analysis']==4167][['TIME_RELATIVE_TO_ALANINE', 'peak_order', 'compound', 'amplitude_2', 'area_2', 'bio_replicate_id', 'TENTATIVE_COMPOUND', 'color']]
    valine_time = int(source_of_truth[source_of_truth['compound'] == 'Val']['TIME_RELATIVE_TO_ALANINE'].values[0])

    possible_valines = df[(df['TIME_RELATIVE_TO_ALANINE'] <= valine_time + TIME_VARIANCE) & (df['TIME_RELATIVE_TO_ALANINE'] >= valine_time - TIME_VARIANCE)]
    df.loc[possible_valines.index, ['TENTATIVE_COMPOUND']] = VAL
    return df

def set_leucines(df):
    source_of_truth = df[df['Analysis']==4167][['TIME_RELATIVE_TO_ALANINE', 'peak_order', 'compound', 'amplitude_2', 'area_2', 'bio_replicate_id', 'TENTATIVE_COMPOUND', 'color']]
    leucine_time = int(source_of_truth[source_of_truth['compound'] == 'Leu']['TIME_RELATIVE_TO_ALANINE'].values[0])
    possible_leucine = df[(df['TIME_RELATIVE_TO_ALANINE'] <= leucine_time + 2) & (df['TIME_RELATIVE_TO_ALANINE'] >= leucine_time - 3)]
    df.loc[possible_leucine.index, ['TENTATIVE_COMPOUND']] = LEU
    return df

def set_isoleucines(df):
    TIME_VARIANCE = 5
    source_of_truth = df[df['Analysis']==4167][['TIME_RELATIVE_TO_ALANINE', 'peak_order', 'compound', 'amplitude_2', 'area_2', 'bio_replicate_id', 'TENTATIVE_COMPOUND', 'color']]
    isoleucine_time = int(source_of_truth[source_of_truth['compound'] == 'Ile']['TIME_RELATIVE_TO_ALANINE'].values[0])
    possible_isoleucine = df[(df['TIME_RELATIVE_TO_ALANINE'] <= isoleucine_time + TIME_VARIANCE) & (df['TIME_RELATIVE_TO_ALANINE'] >= isoleucine_time - TIME_VARIANCE) & (df['amplitude_2'] >= 600)]
    df.loc[possible_isoleucine.index, ['TENTATIVE_COMPOUND']] = ILE
    # this is basically the same functionality as the two lines above -- just taking the max instead of setting a cut off
    # possible_isoleucine = df[(df['TIME_RELATIVE_TO_ALANINE'] <= isoleucine_time + 5) & (df['TIME_RELATIVE_TO_ALANINE'] >= isoleucine_time - 5)]
    # isoleucine_index_ids = []
    # for group, data in possible_isoleucine.groupby('Analysis'):
    #     if len(data) ==1:
    #         isoleucine_index_ids.append(data.index[0])
            
    #     else:
    #         isoleucine_index_ids.append(data['amplitude_2'].idxmax())    
    # df.loc[isoleucine_index_ids, ['TENTATIVE_COMPOUND']] = ILE
    return df

def set_prolines(df):
    TIME_VARIANCE = 5
    source_of_truth = df[df['Analysis']==4167][['TIME_RELATIVE_TO_ALANINE', 'peak_order', 'compound', 'amplitude_2', 'area_2', 'bio_replicate_id', 'TENTATIVE_COMPOUND', 'color']]
    proline_time = int(source_of_truth[source_of_truth['compound'] == 'Pro']['TIME_RELATIVE_TO_ALANINE'].values[0])
    possible_proline = df[(df['TIME_RELATIVE_TO_ALANINE'] <= proline_time + TIME_VARIANCE) & (df['TIME_RELATIVE_TO_ALANINE'] >= proline_time - TIME_VARIANCE)]
    df.loc[possible_proline.index, ['TENTATIVE_COMPOUND']] = PRO
    return df

def set_aspartic_acid(df):
    source_of_truth = df[df['Analysis']==4167][['TIME_RELATIVE_TO_ALANINE', 'peak_order', 'compound', 'amplitude_2', 'area_2', 'bio_replicate_id', 'TENTATIVE_COMPOUND', 'color']]
    aspartic_acid_time = int(source_of_truth[source_of_truth['compound'] == 'Asp']['TIME_RELATIVE_TO_ALANINE'].values[0])
    possible_asp_acid = df[(df['TIME_RELATIVE_TO_ALANINE'] <= aspartic_acid_time + 9) & (df['TIME_RELATIVE_TO_ALANINE'] >= aspartic_acid_time -3) & (df['amplitude_2'] >= 700)]
    df.loc[possible_asp_acid.index, ['TENTATIVE_COMPOUND']] = ASP
    return df

def set_phenylalanine(df):
    source_of_truth = df[df['Analysis']==4167][['TIME_RELATIVE_TO_ALANINE', 'peak_order', 'compound', 'amplitude_2', 'area_2', 'bio_replicate_id', 'TENTATIVE_COMPOUND', 'color']]
    phenylalanine_time = int(source_of_truth[source_of_truth['compound'] == 'Phe']['TIME_RELATIVE_TO_ALANINE'].values[0])
    possible_phenylalanine = df[(df['TIME_RELATIVE_TO_ALANINE'] <= phenylalanine_time + 4) & (df['TIME_RELATIVE_TO_ALANINE'] >= phenylalanine_time - 6) & (df['amplitude_2'] >= 500)]
    df.loc[possible_phenylalanine.index, ['TENTATIVE_COMPOUND']] = PHE
    return df

def set_glutamic_acid(df):
    source_of_truth = df[df['Analysis']==4167][['TIME_RELATIVE_TO_ALANINE', 'peak_order', 'compound', 'amplitude_2', 'area_2', 'bio_replicate_id', 'TENTATIVE_COMPOUND', 'color']]
    glutamic_acid_time = int(source_of_truth[source_of_truth['compound'] == 'Glu']['TIME_RELATIVE_TO_ALANINE'].values[0])
    possible_glutamic_acid = df[(df['TIME_RELATIVE_TO_ALANINE'] <= glutamic_acid_time + 2) & (df['TIME_RELATIVE_TO_ALANINE'] >= glutamic_acid_time -6)]
    df.loc[possible_glutamic_acid.index, ['TENTATIVE_COMPOUND']] = GLU
    return df

def set_lysine(df):
    TIME_VARIANCE = 4
    source_of_truth = df[df['Analysis']==4167][['TIME_RELATIVE_TO_ALANINE', 'peak_order', 'compound', 'amplitude_2', 'area_2', 'bio_replicate_id', 'TENTATIVE_COMPOUND', 'color']]
    lysine_time = int(source_of_truth[source_of_truth['compound'] == 'Lys']['TIME_RELATIVE_TO_ALANINE'].values[0])
    possible_lysine = df[(df['TIME_RELATIVE_TO_ALANINE'] <= lysine_time + TIME_VARIANCE) & (df['TIME_RELATIVE_TO_ALANINE'] >= lysine_time -TIME_VARIANCE)]

    df.loc[possible_lysine.index, ['TENTATIVE_COMPOUND']] = LYS
    return df

def set_tyrosine(df):
    TIME_VARIANCE = 8
    source_of_truth = df[df['Analysis']==4167][['TIME_RELATIVE_TO_ALANINE', 'peak_order', 'compound', 'amplitude_2', 'area_2', 'bio_replicate_id', 'TENTATIVE_COMPOUND', 'color']]
    tyrosine_time = int(source_of_truth[source_of_truth['compound'] == 'Tyr']['TIME_RELATIVE_TO_ALANINE'].values[0])
    possible_tyrosine = df[(df['TIME_RELATIVE_TO_ALANINE'] <= tyrosine_time + TIME_VARIANCE) & (df['TIME_RELATIVE_TO_ALANINE'] >= tyrosine_time -TIME_VARIANCE)]

    df.loc[possible_tyrosine.index, ['TENTATIVE_COMPOUND']] = TYR
    return df


def label_required_aas(df):
    df = set_valines(df)
    df = set_leucines(df)
    df = set_isoleucines(df)
    df = set_prolines(df)
    df = set_aspartic_acid(df)
    df = set_phenylalanine(df)
    return df


def set_methionine(df):
    df = df.groupby('Analysis', group_keys=False).apply(create_time_relative_to_glutamic_acid)
    df = df.groupby('Analysis', group_keys=False).apply(create_time_relative_to_phenylalanine)


def label_non_required_aas(df):
    df = set_glutamic_acid(df)
    df = set_lysine(df)
    df = set_tyrosine(df)
    return df