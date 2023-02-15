import sys

def wt():
    args = sys.argv
    i = args[1]
    while i<args[2]:

        print(i)
        i = i+1
    print("end") 


wt(1,3)

