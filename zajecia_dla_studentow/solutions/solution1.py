characters = "ACDFAAEHHHBBAA?;<<CDEAD98<;;<ADE??<874873012,/."

'''for char in characters:
    out = ord ( str(char) ) - 33
    print(out, end = ' ')
''' 
ph33 = [ord(char) - 33 for char in characters]
ph64 = [ord(char) - 64 for char in characters]
if any(i <0 for i in ph33):
    if all( i > 0 for i in ph64):
        print("Jakość tej sekwencji zapisana jest w systemie phred +64")
        print(ph64)
    else:
        print("Jakość tej sekwencji nie odpowiada ani skali phred +33 ani +64")
else:
    print("Jakość tej sekwencji zapisana jest w systemie phred +33")
    print(ph33)


