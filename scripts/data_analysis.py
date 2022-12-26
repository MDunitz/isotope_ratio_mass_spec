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