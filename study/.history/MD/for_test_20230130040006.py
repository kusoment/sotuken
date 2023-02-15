import sys
import subprocess
import os


args = sys.argv
x = int(args[1])
y = int(args[2])
z = int(args[3])
a = args[4]
V = 0
print(x,y,z,a)

#x = 0
#y = 1
#z = 100

while x<y:
    print(x)
    V = 1000*(x/z*2*3.14159265358979)
    print(V)
    command = ["python.exe","LaplaceCylElFuncFinal.py","3880", "1000", "500", "2000", "0.5", "200", "4", "1000", "20", "1", "0", "Electrode-D5G.csv", V]
    subprocess.run(command)
    
    os.rename("field.npy", str(x)+"-"+str(y) + ".npy")

    x += 1
print("end") 




