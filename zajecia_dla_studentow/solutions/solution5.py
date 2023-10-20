import matplotlib.pyplot as plt

x = [i+1 for i in range(1,len(ph33)+1)]

plt.bar(x,ph33)
plt.ylabel('jakość')
plt.xlabel('nukleotyd')
plt.show()
