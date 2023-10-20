import os
os.system("cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 15 -o cutadapt_out1.fastq data/single-end.fastq")
os.system("cutadapt -a TGGAATTCTCGGGTGCCAAGG -O 8 --trimmed-only -m 15 -o cutadapt_out2.fastq data/single-end.fastq")
# zliczanie liczby rekordow wc -l i dzielenie przez 4

