import sys


args = sys.argv
x = int(args[1])
y = int(args[2])
z = int(args[3])
V = 0


while x<y:
    x = 1
    y = 2
    z = 100
    print(x)
    V = 1000*(x/y*2*3.14159265358979)
    print(V)

    x += 1
print("end") 




