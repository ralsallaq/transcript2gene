#!/usr/bin/env python

import argparse
import os
import sys, time
import numpy as np
import pandas as pd
import glob


def main():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--batches_dir','-d', help='directory where the batches file reside', required=True)
    args_parser.add_argument('--batches_prefix','-p', help='prefix of batch files (e.g. batch_)', required=True)
    args_parser.add_argument('--out_file','-o', help='output file for combined batches', required=True)
    args_parser.add_argument('--failedlist', '-f', help='output failed to pull transcripts', required=True)

    args = args_parser.parse_args()

    batch_dir = args.batches_dir

    prefix = args.batches_prefix

    out_file = args.out_file

    failedlistFile = args.failedlist


    batches = glob.glob(batch_dir+'/'+prefix+'*.tsv')

    batchDFs = [pd.read_csv(f,sep='\t', header=None) for f in batches]

    combined = pd.concat(batchDFs, axis=0)
    
    combined.columns = ['transc','gene','synonyms']

    print(combined.head())

    print(combined.shape)

    combined.to_csv(out_file, sep="\t", header=None, index=False)

    ### remaining transcripts that failed
    failedList = combined[combined['gene'].isnull()]
    
    failedList[['transc']].to_csv(failedlistFile, sep="\t", header=None, index=False)
    

if __name__ == "__main__":
    main()

