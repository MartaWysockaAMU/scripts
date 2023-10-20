#odpowiedzi 2:
import os
import subprocess
import sys
import matplotlib.pyplot as plt

'''
Napisz skrypt w pythonie (moga być 4 skrypty i możesz w nich wykonywać polecenia linuxowe
poprzez np os.system()), który, osobno dla plików wynikowych PE i SE:
2.1. wybierze tylko zmapowane segmenty i tam gdzie występuje hard clipping wykona rozkład
długości przyciętych sekwencji (z końca 5’ i 3’ odczytu) dla każdego mapera; Wykresy zamieść
w protokole;
2.2. wybierze tylko zmapowane segmenty i tam gdzie występuje soft clipping wykona rozkład
długości przyciętych sekwencji (z końca 5’ i 3’ odczytu) dla każdego mapera; Wykresy zamieść
w protokole;

'''

cigar_match = ["M","=","X","I"]

def find_clipping(clipping, cigar):
    idx = cigar.index(clipping)
    s = ""
    for c in range(idx-1,-1,-1):
        if cigar[c].isalpha():
            break
        else:
            if s == "":
                s = cigar[c]
            else:
                s = cigar[c] + s
    return int(s), idx

def count_cigar(cigar):
    x = []
    i = ""
    for c in cigar:
        if c.isalpha():
            if c in cigar_match:
                value = int(i)
                x.append(value)
                i = ""
            else:
                i = ""
        else:
            i += c
    return x

def draw_plot(x, y, xlab, ylab, title):
    plt.bar(x,y)
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.title(title)
    plt.show()

directory = 'OUT/PE/'

b1, b2, bwa, h = [], [], [] , []
#grep ^[^@] b2_1.sam | grep -v -P "SRR[0-9]+\.[0-9]+\t4\t"
for filename in os.listdir(directory):
    bowtie, bowtie2, bwa_bool, hisat2 = False, False, False, False
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        if "b1_" in filename:
            bowtie = True
        elif "b2_" in filename:
            bowtie2 = True
        elif "h_" in filename:
            hisat2 = True
        elif "bwa_" in filename:
            bwa_bool = True
        with open(f) as input:
            for line in input:
                if not line.startswith("@"):
                    if line.strip().split("\t")[2] != "*":
                        if bowtie:
                            b1.append([line.strip().split("\t")[5],line.strip().split("\t")[9]])
                        elif bowtie2:
                            b2.append([line.strip().split("\t")[5],line.strip().split("\t")[9]])
                        elif hisat2:
                            h.append([line.strip().split("\t")[5],line.strip().split("\t")[9]])
                        elif bwa_bool:
                            bwa.append([line.strip().split("\t")[5],line.strip().split("\t")[9]])
                            

b1_mapped_length = {}
b1_S = {}
b1_H = {}
for alignment in b1:
    cigar = alignment[0]
    seq = alignment[1]
    x = count_cigar(cigar)
    if sum(x) in b1_mapped_length:
            b1_mapped_length[sum(x)] += 1
    else:
        b1_mapped_length[sum(x)] = 1
    if "S" in cigar:
        (s, idx_s) = find_clipping("S",cigar)
        if s in b1_S:
            b1_S[s] += 1
        else:
            b1_S[s] =1
        tmp_cigar = cigar[idx_s+1:]
        if "S" in tmp_cigar:
            s, idx_s = find_clipping("S",tmp_cigar)
            if s in b1_S:
                b1_S[s] += 1
            else:
                b1_S[s] =1
    if "H" in cigar:
        (h, idx_h) = find_clipping("H",cigar)
        if h in b1_H:
            b1_H[h] += 1
        else:
            b1_H[h] =1
        tmp_cigar = cigar[idx_h+1:]
        if "S" in tmp_cigar:
            h, idx_h = find_clipping("H",tmp_cigar)
            if h in b1_H:
                b1_H[h] += 1
            else:
                b1_H[h] =1
    
                            
h_mapped_length = {}
h_S = {}
h_H = {}
for alignment in h:
    cigar = alignment[0]
    seq = alignment[1]
    x = count_cigar(cigar)
    if sum(x) in h_mapped_length:
        h_mapped_length[sum(x)] += 1
    else:
        h_mapped_length[sum(x)] = 1
    if "S" in cigar:
        (s, idx_s) = find_clipping("S",cigar)
        if s in h_S:
            h_S[s] += 1
        else:
            h_S[s] =1
        tmp_cigar = cigar[idx_s+1:]
        if "S" in tmp_cigar:
            s, idx_s = find_clipping("S",tmp_cigar)
            if s in h_S:
                h_S[s] += 1
            else:
                h_S[s] =1
    if "H" in cigar:
        (h, idx_h) = find_clipping("H",cigar)
        if h in h_H:
            h_H[h] += 1
        else:
            h_H[h] =1
        tmp_cigar = cigar[idx_h+1:]
        if "S" in tmp_cigar:
            h, idx_h = find_clipping("H",tmp_cigar)
            if h in h_H:
                h_H[h] += 1
            else:
                h_H[h] =1
        
bwa_mapped_length = {}
bwa_S = {}
bwa_H = {}
for alignment in bwa:
    cigar = alignment[0]
    seq = alignment[1]
    x = count_cigar(cigar)    
    if sum(x) in bwa_mapped_length:
        bwa_mapped_length[sum(x)] += 1
    else:
        bwa_mapped_length[sum(x)] = 1
    if "S" in cigar:
        (s, idx_s) = find_clipping("S",cigar)
        if s in bwa_S:
            bwa_S[s] += 1
        else:
            bwa_S[s] =1
        tmp_cigar = cigar[idx_s+1:]
        if "S" in tmp_cigar:
            s, idx_s = find_clipping("S",tmp_cigar)
            if s in bwa_S:
                bwa_S[s] += 1
            else:
                bwa_S[s] =1
    if "H" in cigar:
        (h, idx_h) = find_clipping("H",cigar)
        if h in bwa_H:
            bwa_H[h] += 1
        else:
            bwa_H[h] =1
        tmp_cigar = cigar[idx_h+1:]
        if "S" in tmp_cigar:
            h, idx_h = find_clipping("H",tmp_cigar)
            if h in bwa_H:
                bwa_H[h] += 1
            else:
                bwa_H[h] =1

b2_mapped_length = {}
b2_S = {}
b2_H = {}
for alignment in b2:
    cigar = alignment[0]
    seq = alignment[1]
    x = count_cigar(cigar)
    if sum(x) in b2_mapped_length:
        b2_mapped_length[sum(x)] += 1
    else:
        b2_mapped_length[sum(x)] = 1
    if "S" in cigar:
        (s, idx_s) = find_clipping("S",cigar)
        if s in b2_S:
            b2_S[s] += 1
        else:
            b2_S[s] =1
        tmp_cigar = cigar[idx_s+1:]
        if "S" in tmp_cigar:
            s, idx_s = find_clipping("S",tmp_cigar)
            if s in b2_S:
                b2_S[s] += 1
            else:
                b2_S[s] =1
    if "H" in cigar:
        (h, idx_h) = find_clipping("H",cigar)
        if h in b2_H:
            b2_H[h] += 1
        else:
            b2_H[h] =1
        tmp_cigar = cigar[idx_h+1:]
        if "S" in tmp_cigar:
            h, idx_h = find_clipping("H",tmp_cigar)
            if h in b2_H:
                b2_H[h] += 1
            else:
                b2_H[h] =1
        
x = list(h_mapped_length.keys())
y = [h_mapped_length[k] for k in x]
draw_plot(x, y, 'Długość zmapowania', 'Liczba zmapowań', 'Hisat2 - alignment')

x = list(b1_mapped_length.keys())
y = [b1_mapped_length[k] for k in x]
draw_plot(x, y, 'Długość zmapowania', 'Liczba zmapowań', 'Bowtie - alignment')

x = list(b2_mapped_length.keys())
y = [b2_mapped_length[k] for k in x]
draw_plot(x, y, 'Długość zmapowania', 'Liczba zmapowań', 'Bowtie2 - alignment')

x = list(bwa_mapped_length.keys())
y = [bwa_mapped_length[k] for k in x]
draw_plot(x, y, 'Długość zmapowania', 'Liczba zmapowań', 'BWA - alignment')

x = list(h_S.keys())
y = [h_S[k] for k in x]
draw_plot(x, y, 'Długość soft clippingu', 'Liczba odczytów', 'Hisat2 - soft clipping')

x = list(b1_S.keys())
y = [b1_S[k] for k in x]
draw_plot(x, y, 'Długość soft clippingu', 'Liczba odczytów', 'Bowtie - soft clipping')

x = list(b2_S.keys())
y = [b2_S[k] for k in x]
draw_plot(x, y, 'Długość soft clippingu', 'Liczba odczytów', 'Bowtie2 - soft clipping')

x = list(bwa_S.keys())
y = [bwa_S[k] for k in x]
draw_plot(x, y, 'Długość soft clippingu', 'Liczba odczytów', 'Bwa - soft clipping')

x = list(h_H.keys())
y = [h_H[k] for k in x]
draw_plot(x, y, 'Długość soft clippingu', 'Liczba odczytów', 'Hisat2 - hard clipping')

x = list(b1_H.keys())
y = [b1_H[k] for k in x]
draw_plot(x, y, 'Długość soft clippingu', 'Liczba odczytów', 'Bowtie - hard clipping')

x = list(b2_H.keys())
y = [b2_H[k] for k in x]
draw_plot(x, y, 'Długość soft clippingu', 'Liczba odczytów', 'Bowtie2 - hard clipping')

x = list(bwa_H.keys())
y = [bwa_H[k] for k in x]
draw_plot(x, y, 'Długość soft clippingu', 'Liczba odczytów', 'Bwa - hard clipping')

