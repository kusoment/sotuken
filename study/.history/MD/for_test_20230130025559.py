import sys

args = sys.argv
x = int(args[1])
y = int(args[2])
z = int(args[3])
V = 0


while x<y:

    print(x)
    V = 1000*(x/y)

    x += 1
print("end") 




