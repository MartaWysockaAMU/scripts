import numpy as np

out_list = []
for char in characters:
    out = ord ( str(char) ) - 33
    out_list.append(out)

mini = min(out_list)
maxi = max(out_list)
'''
if sum([i == maxi for i in ph33]) > 1:
    print("Nukleotydy o najwyższej jakosci znajdują sie na pozycjach {}".format([i + 1 for i in np.where(np.array(ph33) == maxi)[0]]))
else:
    print("Nukleotyd o najwyższej jakosci znajduje sie na pozycji {}".format( out_list.index( maxi)+1 ))

if sum([i == mini for i in ph33]) > 1:
    print("Nukleotydy o najniższej jakosci znajdują sie na pozycjach {}".format([i + 1 for i in np.where(np.array(ph33) == mini)[0]]))
else:
    print("Nukleotyd o najnizszej jakosci znajduje sie na pozycji {}".format( out_list.index( mini)+1 ))
'''

print("Nukleotyd o najnizszej jakosci znajduje sie na pozycji {}".format( out_list.index( mini)+1 ))
print("Nukleotyd o najwyższej jakosci znajduje sie na pozycji {}".format( out_list.index( maxi)+1 ))

