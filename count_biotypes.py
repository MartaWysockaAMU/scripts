import os
import sys
import pysam

mRNA = [] # chromosom, pozycje start i stop z pliku gff dla genów gdzie biotype=protein_coding 
#grep biotype=protein_coding ../Ribo_benchmark_data/genome_saureus/GCF_000013425.1_ASM1342v1_genomic.gff
rRNA = [] # chromosom, pozycje start i stop z pliku gff dla genów gdzie w 3ciej kolumnie jest to rRNA oraz ID=rna
# grep -P '\trRNA\t' ../Ribo_benchmark_data/genome_saureus/GCF_000013425.1_ASM1342v1_genomic.gff
tRNA = [] # chromosom, pozycje start i stop z pliku gff dla genów gdzie product=tRNA a 3kolumna to exon - zawierają się wtedy tRNA i pseudogenic tRNA
#grep tRNA ../Ribo_benchmark_data/genome_saureus/GCF_000013425.1_ASM1342v1_genomic.gff | grep -P '\texon\t'
done = []

cmd = "grep biotype=protein_coding ../Ribo_benchmark_data/genome_saureus/GCF_000013425.1_ASM1342v1_genomic.gff > mRNA.gff"
os.system(cmd)
with open("mRNA.gff") as gfffile:
    for line in gfffile:
        tmp = line.strip().split('\t')
        coordinants = [tmp[0],int(tmp[3]),int(tmp[4])]
        if coordinants not in mRNA:
            mRNA.append(coordinants)

cmd = "grep -P '\trRNA\t' ../Ribo_benchmark_data/genome_saureus/GCF_000013425.1_ASM1342v1_genomic.gff > rRNA.gff"
os.system(cmd)
with open("rRNA.gff") as gfffile:
    for line in gfffile:
        tmp = line.strip().split('\t')
        coordinants = [tmp[0],int(tmp[3]),int(tmp[4])]
        if coordinants not in mRNA:
            mRNA.append(coordinants)

cmd = "grep tRNA ../Ribo_benchmark_data/genome_saureus/GCF_000013425.1_ASM1342v1_genomic.gff | grep -P '\texon\t' > tRNA.gff"
os.system(cmd)
with open("tRNA.gff") as gfffile:
    for line in gfffile:
        tmp = line.strip().split('\t')
        coordinants = [tmp[0],int(tmp[3]),int(tmp[4])]
        if coordinants not in mRNA:
            mRNA.append(coordinants)

c=0
d = 'mapping/STAR/'
for filename in os.listdir(d):
    if '.bam' in filename and not '.bai' in filename:
        f = os.path.join(d, filename)
        print(f)
        #with open(f) as file:
        mybam = pysam.AlignmentFile(f)
        for m in rRNA:
            iter = mybam.fetch(m[0],m[1],m[2])
            for x in iter:
                
                if (x.tags[0][1]) > 1:
                    print(str(x))
                c += 1
                #print(x.query_name)
                if x.query_name not in done:
                    done.append(x.query_name)
                    sys.exit()
            print(len(done), c)
            sys.exit()

