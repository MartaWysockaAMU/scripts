strings = {} #seq
prevLine = ""
k=3
with open("rosalind_grph.txt") as input_file:
    seq = ""
    prevHeader = ""
    for line in input_file:
        if line.startswith(">"):
            if prevHeader == "":
                prevHeader = line.strip()[1:]
            else:
                suf = seq[-k:]
                pref = seq[0:k]
                strings[prevHeader] = [seq,pref,suf]
            prevHeader = line.strip()[1:]
            seq = ""
        else:
            seq += line.strip()
    suf = seq[-k:]
    pref = seq[0:k]
    strings[prevHeader] = [seq,pref,suf]
            
graph = []
for s in strings:
    head_s = strings[s][1]
    tail_s = strings[s][2]
    for t in strings:
        head_t = strings[t][1]
        tail_t = strings[t][2] 
        if s != t:
            edge = None
            if tail_s == head_t: # s->t
                edge = s+" "+t
            elif tail_t == head_s: #t->s
                edge = t+" "+s
            if edge and (edge not in graph):
                    graph.append(edge)

for g in graph:
    print(g)