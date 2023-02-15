import sys


args = sys.argv
#x = int(args[1])
#y = int(args[2])
#z = int(args[3])
V = 0

x = 0
y = 51
z = 100

while x<y:
    print(x)
    V = 1000*(x/y*2*3.14159265358979)
    V2 = ('{:.5g}'.format(V))
    print(V2)

    x += 1
print("end") 




