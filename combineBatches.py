#!/usr/bin/env python
import argparse
import os
import sys, time
import numpy as np
import pandas as pd
import glob


def main():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--accessions_file', '-a', help='the original accession file from vep cache', required=True) 
    args_parser.add_argument('--batches_dir','-d', help='directory where the batches file reside', required=True)
    args_parser.add_argument('--batches_prefix','-p', help='prefix of batch files (e.g. batch_)', required=True)
    args_parser.add_argument('--out_file','-o', help='output file for combined batches', required=True)
    args_parser.add_argument('--failedlist', '-f', help='output failed to pull transcripts', required=True)

    args = args_parser.parse_args()

    original_accessions_file = args.accessions_file

    batch_dir = args.batches_dir

    prefix = args.batches_prefix

    out_file = args.out_file

    failedlistFile = args.failedlist

    orig_accession = pd.read_csv(original_accessions_file, sep="\t", header=None)
    orig_accession.columns = ['accession'] 
    orig_accession.loc[:,'base'] = orig_accession['accession'].apply(lambda r: r.split(".")[0]) 
 
    batches = glob.glob(batch_dir+'/'+prefix+'*.tsv')

    batchDFs = [pd.read_csv(f,sep='\t', header=None) for f in batches]

    combined = pd.concat(batchDFs, axis=0)
    
    combined.columns = ['transc','gene','synonyms']

    print(combined.head())

    print(combined.shape)


    ### remaining transcripts that failed
    failedList = combined[combined['gene'].isnull()]
    failedList.loc[:,'base'] = failedList['transc'].apply(lambda r: r.split(".")[0].strip()) 
    failedList = failedList[['base']]

    tr_idx = orig_accession['base'].isin(failedList['base'].values)
    tr = orig_accession.loc[tr_idx].set_index('base')['accession'].to_dict()
    failedList = failedList.replace(tr)

    combined = combined[~combined['gene'].isnull()]

    combined.to_csv(out_file, sep="\t", header=None, index=False)
    failedList.to_csv(failedlistFile, sep="\t", header=None, index=False)
    

if __name__ == "__main__":
    main()

