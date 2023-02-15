import sys

args = sys.argv
i = args[1]
print(i)
while i<args[2]:
    print(i)
    i += 1

print("end")