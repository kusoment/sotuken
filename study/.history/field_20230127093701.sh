#!/bin/bash

i=$1
V=0
while [ $i -lt $2 ] #第一引数　<　第二引数
do
    echo $i
    V=`echo "scale=5; 1000 * c (  $i / $3 * 2 * 3.14159265358979 )" | bc -l` #scale 小数位 
    mpirun -n 4 python ./LaplaceCylElFuncFinal.py 3880 1000 500 2000 0.5 200 4 1000 20 1 0 $4 $V 
    mv field.npy "$i-$3.npy" #rename $-100 ($+1)-100...と作られる($2まで)
    i=`expr $i + 1`
done