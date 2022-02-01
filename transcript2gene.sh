#first load python
module load python/3.7.0
nthreads=100
mkdir logging
#bsub -n 1 -R "span[hosts=1] rusage[mem=500]" -P preferredIso -J "ncbi[1-$nthreads]" -o "logging/output.%I" "python get_geneSymbolsFromAccessions.py --batches -i refseqs_100_GRCh38.txt -e ralsalla@stjude.org  -t '$nthreads' -r \$LSB_JOBINDEX"
bsub -n 1 -R "span[hosts=1] rusage[mem=500]" -P preferredIso -J "ncbi[1-$nthreads]" -o "logging/output.%I" "python get_geneSymbolsFromAccessions.py --batches --update -i vep_transcripts_genes_catchup.tsv -e ralsalla@stjude.org  -t '$nthreads' -r \$LSB_JOBINDEX"
