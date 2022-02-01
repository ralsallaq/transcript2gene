#!/usr/bin/env python
import argparse
import os
import sys
import time
import numpy as np
import pandas as pd
from Bio import Entrez


def get_id(accession):
    """ seaches ncbi refseq Gene database for the transcript and returns id """
    idfound = False
    term = 'srcdb_refseq[property] AND "Official Symbol" AND '+str(accession)
    while not idfound:
        try:

            time.sleep(8)
            search = Entrez.esearch(term=term, db="Gene", retmode="xml")
            record = Entrez.read(search)
            idfound = True

        except:

            idfound = False
            #print("searching again for ", term, " in 8''", file=sys.stderr)
            time.sleep(8)


    if len(record['IdList']) > 0:
        return record['IdList'][0]
    else:
        return None

def fetch_data(iid):
    """ uses an id to get data """
    data_fetched = False
    while not data_fetched:
        try:
            time.sleep(8)
            fetch = Entrez.efetch(id=iid, db="Gene", retmode="xml")
            data = Entrez.read(fetch)
            data_fetched = True
        except:
            data_fetched = False
            #print("trying again for iid", iid, " in 8''", file=sys.stderr)
            time.sleep(8)
    official_symbol = data[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']

    try:

        synonyms = data[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']

    except KeyError:

        synonyms = None

    return official_symbol, synonyms


def main():
    """
       Given accession(s) for transcripts this script finds the official symbol of the gene
       and all synonyms (including official symbol)
       Usage:
       get_geneSymbolsFromAccessions.py \
            -i <list_of_accessions_file> \
            -o <accessions_and_gene_symboles_file>

    """
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--accessions_file', '-i', 
            help='headerless file that encompasses a vertical list of accessions for transcripts', required=True)
    args_parser.add_argument('--yourEmail', '-e', 
            help='a valid email for NCBI server inquiries', required=True)
    args_parser.add_argument('--out_file', '-o', 
            help='output file in TSV format', required=not '--batches' in sys.argv)
    args_parser.add_argument('--batches', 
            help='run in batches; can be used when the number of accession is large, defaule=False', action='store_true')
    args_parser.add_argument('--num_threads', '-t', 
            type=int, help='number of threads; required with --batches', required ='--batches' in sys.argv)
    args_parser.add_argument('--batch_r', '-r', type=int, 
            help='the batch run number (1,2, ..., num_threads); required with --batches and must be <= num_threads', required='--batches' in sys.argv )
    args_parser.add_argument('--update', 
            help='update transcripts to a newer version before pulling info from NCBI; defaule=False', action='store_true')

    args = args_parser.parse_args()

    # First handle the files
    accesssions_file = args.accessions_file

    ### set up the email for NCBI inquiries 
    Entrez.email = args.yourEmail

    #print("............read accessions file............", file=sys.stderr)
    ##### read the list of accessions
    access = pd.read_csv(accesssions_file, header=None)

    if args.batches:

        batchN = args.batch_r 

        num_of_threads = args.num_threads 

        assert((batchN >= 1) and (batchN <= num_of_threads)), "invalid batch number! batch number must be a positive integer (1,2,3,...NumThreads)"

        #temp_dir = Path(tempfile.TemporaryDirectory().name)

        #print("creating a batches directory", file=sys.stderr)

        os.makedirs("./batches", exist_ok=True)


        #split the dataframe into a list of dataframes (batches)
        list_df_batches = np.array_split(access, num_of_threads)

        batch_df = list_df_batches[batchN-1]

        #print("batch Number", batchN, file=sys.stderr)
        #print("batch DF head", batch_df.head(), file=sys.stderr)

        output_file = './batches/batch_df_'+str(batchN)+'.tsv'
        #print(output_file, file=sys.stderr)

    else:

        batch_df = access
        output_file = args.out_file


    assert(batch_df.columns.isin([0]).sum() == 1), "the accessions file must be headerless"

    #print(batch_df.head(), file=sys.stderr)

    ### set up info to None
    batch_df.loc[:, 'Gene'] = None
    batch_df.loc[:, 'Synonym'] = None

    for i, acc in batch_df.iterrows():

        if args.update:

            updated_acc = acc[0]

            num_of_tries = 0

            iid = None

            while (iid is None)  and (num_of_tries < 5):

                updated_acc = updated_acc.split(".")[0]+'.'+str(int(acc[0].split(".")[-1])+1)

                #print("the accession ", acc[0], " is updated to ", updated_acc, " before lookup", file=sys.stderr)

                batch_df.loc[i, 0] = updated_acc

                iid = get_id(updated_acc)

                num_of_tries += 1


        else:

            iid = get_id(acc[0])

        if iid:

            official_symbol, synonyms = fetch_data(iid)

            #print(i, acc[0], official_symbol, synonyms, file=sys.stderr)

            batch_df.loc[i, 'Gene'] = official_symbol 

            batch_df.loc[i, 'Synonym'] = ";".join(synonyms) if synonyms else None


    batch_df.columns = ['accession', 'gene', 'synonyms']
    batch_df.to_csv(output_file, sep="\t", index=False, header=None)
    
    
if __name__ == "__main__":
    main()
