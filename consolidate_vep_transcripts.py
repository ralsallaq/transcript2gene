# coding: utf-8
import sys
import numpy as np
import pandas as pd

run1 = pd.read_csv("combined_Gene_update.tsv", sep="\t", header=None)
assert(run1[1].isnull().sum() == 0),"missing genes in the Gene's run"
run1 = run1[[0,1]]
run1.columns = ['accession_versioned','gene']

run2 = pd.read_csv("combined_Nucleotide.tsv", sep="\t", header=None)
assert(run2[1].isnull().sum() == 0), "missing genes in the Nucleotide's run"
run2 = run2[[0,1]]
run2.columns = ['accession_versioned','gene']

manual = pd.read_csv("combined_manual.tsv", sep="\t", header=None)
assert(manual[1].isnull().sum() == 0), "missing genes in the maunally obtained"
manual = manual[[0,1]]
manual.columns = ['accession_versioned','gene']

# getting accessions
run1.loc[:,'accession'] = run1['accession_versioned'].apply(lambda r:r.split(".")[0])
run2.loc[:,'accession'] = run2['accession_versioned'].apply(lambda r:r.split(".")[0])
manual.loc[:,'accession'] = manual['accession_versioned'].apply(lambda r:r.split(".")[0])

#### if want to append the updated transcripts 
found_ = run1.append(run2)
found_ = found_.append(manual)
print(found_.shape)
print(found_.tail())

found_ = found_[['gene','accession','accession_versioned']]

found_.to_csv("consolidated_vep_transcripts_genes.tsv", sep="\t", index=False)
