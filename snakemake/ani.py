'''
Title: ani.py
Author: Renee Oles
Date: 7/8/2020
Program: takes text file with bin comparisons based off percent identity, calculates top three bin that are closest in
percent identity for each bin
Input: ani.py -f <file.txt> -o <name of output file>
     optional flags: -p <percent identity cutoff for filtering> -m <minimum number of matches for filtering>
Output: csv file with top three bin matches for each bin
'''
import os
import pandas as pd
import argparse
import numpy as np
import csv


def count_ani(df, output):
    global bin_df
    top = pd.DataFrame()
    for sample in df.Sample1.unique():
        for bin1 in df[df['Sample1'] == sample].Bin1.unique():
            i = 0
            mask = df[(df['Bin1'] == bin1) & (df['Sample1'] == sample)]
            mask = mask[mask['Sample2'] != sample]
            for sample2 in mask.Sample2.unique():
                    for bin2 in mask[mask['Sample2'] == sample2].Bin2.unique():
                        mask2 = mask[(mask['Sample2'] == sample2) & (mask['Bin2'] == bin2)]
                        total = calculate_total(mask2)
                        matches = add_matches(mask2)
                        percent_identity = avg_percent_identity(mask2)
                        if percent_identity >= 80.0:
                            bin_df = pd.DataFrame(
                                {'Sample1': sample, 'Sample2': sample2, 'Bin1': bin1, 'Bin2': bin2, 'Score': total,
                                 'Matches': matches, 'Average Percent Identity': percent_identity}, index=[0])
            #bin_df = bin_df.sort_values(by='Score', ascending=False)
                        top = pd.concat([top, bin_df])
        df = df[df['Sample2'] != sample]
        df = df[df['Sample1'] != sample]
    top.to_csv("%s.csv" % output, index=False)


def calculate_total(df):
    total = 0
    for ind in df.index:
        total = total + df.loc[ind]['Percent_Identity']/100 * df.loc[ind]['Hits']
    return total


def add_matches(df):
    matches = 0
    for ind in df.index:
        matches = matches + df.loc[ind]['Hits']
    return matches


def avg_percent_identity(df):
    percent = 0
    total = 0
    for ind in df.index:
        percent = percent + df.loc[ind]['Percent_Identity']
        total = total + 1
    if total > 0:
        percent = percent/total
    return percent

def trim(df, matches, percentage):
    df = df[(df['Percent_Identity'] >= percentage) & (df['Hits'] >= matches) & (df['Bin1'] != "NoBin") & (df['Bin2'] != "NoBin")]
    return df

def cluster(df, output):
    df['Sample1'] = df['Sample1'] + "_" + df['Bin1']
    df['Sample2'] = df['Sample2'] + "_" + df['Bin2']
    df = df[['Sample1', 'Sample2']]
    list = []
    for i in df.Sample1.unique():
        df_sub = df[df['Sample1'] == i]
        if len(df_sub) > 0:
            line1 = df_sub['Sample2'].str.cat(sep=',')
            line1 = i + "," + line1
            list.append(line1)
        for j in df_sub.Sample2:
            df = df[df['Sample1'] != j]
            df = df[df['Sample2'] != j]
    pd.DataFrame(list).to_csv("%s_cluster.csv" % output, index=None, header=None)


if __name__ == "__main__":
    """ Arguments """
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-f", "--File", help="Txt file",
                        required=True)
    parser.add_argument("-o", "--Output", help="Csv output filename",
                        required=True)
    parser.add_argument("-p", "--Percentage", help="List percentage cutoff",
                        required=False)
    parser.add_argument("-m", "--Matches", help="List number of matches cutoff",
                        required=False)
    argument = parser.parse_args()
    df = pd.read_table(argument.File, dtype={"Percent_Identity": float, "Hits": float})
    if argument.Percentage and argument.Matches:
        df = trim(df, float(argument.Matches), float(argument.Percentage))
    elif argument.Pertage:
        df = trim(df, 0, float(argument.Percentage))
    elif argument.Matches:
        df = trim(df, float(argument.Matches), 0)
    #count_ani(df, argument.Output)
    df = pd.read_csv("%s.csv" % argument.Output)
    df = df[df['Average Percent Identity'] >= 90]
    cluster(df, argument.Output)