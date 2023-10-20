RNA = ""
introns = []
with open("rosalind_splc.txt") as input_file:
    header = next(input_file)
    line = next(input_file)
    #RNA += line.strip()
    while not line.startswith(">"):
        RNA += line.strip()
        line = next(input_file)
    intron = ""
    for line in input_file:
        if line.startswith(">"):
            if intron != "":
                introns.append(intron)
                intron = ""
        else:
            intron += line.strip()
    if intron != "":
        introns.append(intron)
for i in introns:
    RNA = RNA.replace(i,"")

codons = {}
stop = [] # TAA, TAG, TGA
with open("codon_table.txt") as file1:
    for line in file1:
        nt = line.strip().split()[0]
        aa = line.strip().split()[1]
        if aa == "*":
            stop.append(nt)
            codons[nt] = aa
        else:
            codons[nt] = aa
aa = ""
for i in range(0,len(RNA),3):
    codon =""
    for j in range(i,i+3):
        codon += RNA[j]
    if codon in stop:
        break
    else:
        aa += codons[codon]

print(aa)
