#conda activate qc
#whereis fastqc
#ls -l /home/linux/anaconda3/envs/qc/bin/ | grep fastqc
#ls -l /home/linux/anaconda3/envs/qc/opt/fastqc-0.12.1/Configuration/


import os
import statistics as st

file = "data/single-end2.fastq"
adapter_list = "/etc/fastqc/Configuration/adapter_list.txt"
#adapter_list = "/home/linux/anaconda3/envs/qc/opt/fastqc-0.12.1/Configuration/adapter_list.txt"

#os.system("fastqc --extract {}".format(file))
adapters = {}
with open ("data/single-end2_fastqc/fastqc_data.txt") as data:
    for line in data:
        if line.startswith(">>Adapter"):
            header = next(data).strip().split("\t")
            for h in header:
                adapters[h] = []
            line = next(data)
            while not line.startswith(">>"):
                values = line.strip().split("\t")
                for i in range(len(header)):
                    adapters[header[i]].append(float(values[i]))
                line = next(data)
                
adapter_delete = ""
for adapter in adapters:
    m = 0.0
    if adapter != '#Position':
        indexes = [ n for n,i in enumerate(adapters[adapter]) if i>0.0 ]
        if len(indexes) > 0:
            ix = indexes[0]
            tmp = st.mean(adapters[adapter][ix:])
            if tmp > m:
                m = tmp
                adapter_delete = adapter
                
print(adapter_delete)

with open (adapter_list) as l:
    for line in l:
        if adapter_delete in line:
            seq = line.strip().split()[-1]

print("Usuniety powinien zostać adapter {} o sekwencji {}\n".format(adapter_delete, seq))
out_file = "cutadapt_out_7.fastq"
#os.system("cutadapt -a {} -m 15 -o {} data/single-end2.fastq".format(seq, out_file))
#os.system("fastqc {}".format(out_file))
