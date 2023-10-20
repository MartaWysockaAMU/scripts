os.system("fastq_quality_trimmer -t 20 -i data/paired-end_1.fastq -o fastx_out1pe1.fastq")
os.system("fastq_quality_trimmer -t 20 -i data/paired-end_2.fastq -o fastx_out1pe2.fastq")
os.system("fastq_quality_filter -p 90 -q 20 -i fastx_out1pe1.fastq -o fastx_out2pe1.fastq")
os.system("fastq_quality_filter -p 90 -q 20 -i fastx_out1pe2.fastq -o fastx_out2pe2.fastq")
