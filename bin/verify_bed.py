"""
Ugly code to try validate bed file created from nirvana_refseq gff file with gff_to_bed.py

Have to use .gtf from UCSC which has 2 issues:
    - 1 based, nirvana bed is 0 based
    - doesn't include stop codons where nirvana does

This is fixed in gtf_adjust before identifying how many match between the UCSC transcripts and Nirvana


Jethro Rainford 200512
"""

import pandas as pd
import os
import sys

from itertools import tee

def input_files():
    """
    Set input files for UCSC and Nirvana
    Extract just the CDS lines in .gtf downloaded from UCSC with:
    $ awk '$3== "CDS"' ucsc_transcripts.gtf  > ucsc_CDS.gtf

    Args: None

    Returns: 
        ucsc_bed (str): file path to UCSC .gtf file
        nirvana_bed (str): file path to Nirvana bed file
    """

    ucsc_bed = os.path.expanduser("./ucsc_CDS.gtf")
    nirvana_bed = os.path.expanduser("./refseq_nirvana_2.0.10.bed")

    return ucsc_bed, nirvana_bed


def bed_to_df(ucsc_bed, nirvana_bed):
    """
    Make df's of both input files

    Args:
        ucsc_bed (str): file path to UCSC .gtf file
        nirvana_bed (str): file path to Nirvana bed file

    Returns:
        ucsc_df (df): dataframe of UCSC gtf file
        nirvana_df (df): dataframe of Nirvana bed file
    """

    # read in nirvana bed and set header
    nirvana_df = pd.read_csv(nirvana_bed, sep='\t', names=["chrom", "start", "end", "transcript"], low_memory=False)

    # UCSC gtf to df, have to set header later as has lots of columns to remove
    ucsc_df = pd.read_csv(ucsc_bed, sep='\t', header=None)
    
    # get just the transcript name from last column as it's in a long string 
    transcript_col = ucsc_df.iloc[:,-1].str.split('"',0).str[3]
    ucsc_df.iloc[:,-1] = transcript_col

    return ucsc_df, nirvana_df


def gtf_adjust(ucsc_df):
    """
    Function to adjust UCSC gtf for being 1 based and not inc. stop codons

    Args:
        ucsc_df (df): dataframe of UCSC gtf file

    Returns: 
        ucsc_df (df): dataframe of UCSC gtf file with adjusted bases and 
        only required columns for comparison
    """
    # pair iter function
    def pairwise(iterable):
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)
    
    def minus_1(num):
        return int(num) - 1
    
    def strip_chr(chrom):
        return chrom.strip("chr")

    # keep only rows with chr 1-22/X/Y (i.e. remove those with _fix_ etc.)
    print("Checking chromosome entries")

    rows_to_drop = []

    for i, row in ucsc_df.iterrows():
        if "_" in row[0]:
            rows_to_drop.append(i)

    ucsc_df = ucsc_df.drop(ucsc_df.index[rows_to_drop])

    # drop 'chr' prefix from chromosome to match nirvana
    ucsc_df[0] = ucsc_df[0].apply(strip_chr)

    # minus 1 from start of each as gff is 1 based and nirvana is 0 based
    ucsc_df[3] = ucsc_df[3].apply(minus_1)

    print("Adjusting for stop codons")

    for (i1, row1), (i2, row2) in pairwise(ucsc_df.iterrows()):
        # loop to adjust by 3bp of first exon if + strand and last exon if - strand

        transcript = row1[8]
        next_transcript = row2[8]

        strand = row1[6]
        next_strand = row2[6]

        if strand == "+":
            if transcript != next_transcript:
                # if next transcript is diff. => last exon => adjust last by +3
                ucsc_df.iloc[i1, 4] = int(row1[4]) + 3
            else:
                continue

        if next_strand == "-":
            if next_transcript != transcript:
                # next transcript is different & strand of next is -ve => 
                # adjust start of the next transcript by -3
                ucsc_df.iloc[i2, 3] = int(row2[3]) - 3
            else:
                continue
    
    # create output .gff with adjusted bases
    ucsc_df.to_csv("./ucsc_CDS_base_adjusted.gtf", index = None, sep='\t', mode='a')

    # just leave start, end, transcript columns & set header
    ucsc_df = ucsc_df.drop(ucsc_df.columns[[1,2,5,6,7]], axis = 1)
    header = ["chrom", "start", "end", "transcript"]
    ucsc_df.columns = header

    print(ucsc_df)

    return ucsc_df


def calc_diff(nirvana_df, ucsc_df):
    """
    Get differences between UCSC and Nirvana, and write to output files

    Args:
        ucsc_df (df): dataframe of UCSC gtf file
        nirvana_df (df): dataframe of Nirvana bed file
    
    Returns: None
    """

    # find differences between ucsc and nirvana
    diff = pd.merge(ucsc_df, nirvana_df, how="outer",indicator=True)

    only_in_ucsc = diff[diff['_merge'] == "left_only"] # entries only in ucsc df
    only_in_nirvana = diff[diff['_merge'] == "right_only"] # entries only in the nirvana df
    both = diff[diff['_merge'] == "both"] # entries matching between both

    print("only in ucsc")
    print(only_in_ucsc)

    print("in both")
    print(both)

    # write output files
    only_in_ucsc.to_csv("./only_in_ucsc.bed", index=None, sep='\t', mode='a')
    both.to_csv("./ucsc_nirvana_match_transcripts.bed", index=None, sep='\t', mode='a')


if __name__ == "__main__":

    ucsc_bed, nirvana_bed = input_files()

    nirvana_df, ucsc_df = bed_to_df(ucsc_bed, nirvana_bed)

    ucsc_df = gtf_adjust(ucsc_df)

    calc_diff(nirvana_df, ucsc_df)