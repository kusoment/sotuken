import sys
import subprocess
from subprocess import PIPE
import os


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
    command = ["python.exe","LaplaceCylElFuncFinal.py","3880", "1000", "500", "2000", "0.5", "200", "4", "1000", "20", "1", "0", "Electrode-D5G.csv", V2]
    subprocess.run(command)
    
    os.rename("field.npy", str(x)+"-"+str(y) + ".npy")

    x += 1
print("end") 




