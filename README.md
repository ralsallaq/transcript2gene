# transcript2gene
pulls gene symbols for a list of transcripts from NCBI 

```
usage: transcript2gene.py [-h] --accessions_file ACCESSIONS_FILE --yourEmail

                          YOUREMAIL --out_file OUT_FILE [--batches]
                          
                          [--num_threads NUM_THREADS] [--batch_r BATCH_R]
                          
                          [--update]

optional arguments:

  -h, --help            show this help message and exit
  
  --accessions_file ACCESSIONS_FILE, -i ACCESSIONS_FILE
  
                        headerless file that encompasses a vertical list of
                        
                        accessions for transcripts
                        
  --yourEmail YOUREMAIL, -e YOUREMAIL
  
                        a valid email for NCBI server inquiries
                        
  --out_file OUT_FILE, -o OUT_FILE
  
                        output file in TSV format
                        
  --batches             run in batches; can be used when the number of
  
                        accession is large, defaule=False
                        
  --num_threads NUM_THREADS, -t NUM_THREADS
  
                        number of threads; required with --batches
                        
  --batch_r BATCH_R, -r BATCH_R
  
                        the batch run number (1,2, ..., num_threads); required
                        
                        with --batches and must be <= num_threads
                        
  --update              update transcripts to a newer version before pulling
  
                        info from NCBI; defaule=False
                        
```
