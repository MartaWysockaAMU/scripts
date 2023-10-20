from itertools import combinations 
strings = []
seq = ""
with open("rosalind_lcsm.txt") as file1:
    for line in file1:
        if line.startswith(">"):
            if seq != "":
                strings.append(seq)
                seq = ""
        else:
            seq += line.strip()
    if seq != "":
        strings.append(seq)

sorted_strings = sorted(strings, key=len)
all_substrings = []
shortest_length = len(sorted_strings[0])
shortest_string = sorted_strings[0]

all_substrings = [shortest_string[x:y] for x, y in combinations( 
            range(len(shortest_string) + 1), r = 2)] 


all_substrings = sorted(all_substrings, key=len, reverse=True)
#print(all_substrings)
#print(all_substrings[0])
longest_shared = None
for ss in all_substrings:
    candidate = True
    for s in strings:
        if not ss in s:
            candidate = False
            break
    if candidate:
        longest_shared = ss
        break

print(longest_shared)
