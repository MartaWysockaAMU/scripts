os.system("fastq_quality_trimmer -t 20 -l 20 -i cutadapt_out2.fastq -o fastx_out1.fastq")
os.system("fastq_quality_filter -p 90 -q 20 -i fastx_out1.fastq -o fastx_out2.fastq")
