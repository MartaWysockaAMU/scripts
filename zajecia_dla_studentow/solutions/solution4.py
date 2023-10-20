'''out_array = np.asarray(out_list)
print(out_array[out_array<20].shape[0]) # liczba elementow o jakosci ponizej 20
print((out_array[out_array<20].shape[0]/out_array.shape[0])*100) # jaka czesc wszystkich stanowia (w %)'''

print("Jest {} pozycji o jakości poniżej 20 w skali phred".format(len(np.where(np.array(ph33) < 20)[0])))
print("stanowią one {}% wszystkich nukleotydów".format(round(len(np.where(np.array(ph33) < 20)[0])/len(ph33) * 100,2)))
