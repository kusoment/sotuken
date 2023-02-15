
from django.shortcuts import render
from django.http import HttpResponse
from django.views import View
import csv
import pandas as pd
import numpy as np
import subprocess
import os

# Create your views here.
def index(request):
    context = {
        'foo': 'foo string!',
    }
    return render(request, 'MD/index.html', context)


#条件入力
class efView(View):
    def get(self, request, *args, **kwargs):
        context = {
            'message': "電極条件を入力してください",
        }
        return render(request, 'MD/ef.html', context)
    
    def post(self, request, *args, **kwargs):
        context = {
            'message': "送信されました",
        }
        print(request.POST['num1_vol'],request.POST['num2_vol'])
        print(request.POST['num1_vol'])
        #貰った値をpandas,numpyでcsv化　参考https://www.self-study-blog.com/dokugaku/python-pandas-dataframe-make/
        df1 = pd.DataFrame([
            [request.POST['num1_vol'],request.POST['num1_X1'],request.POST['num1_Y1'],request.POST['num1_X2'],request.POST['num1_Y2'],request.POST['num1_RF']],
            [request.POST['num2_vol'],request.POST['num2_X1'],request.POST['num2_Y1'],request.POST['num2_X2'],request.POST['num2_Y2'],request.POST['num2_RF']],
            [request.POST['num3_vol'],request.POST['num3_X1'],request.POST['num3_Y1'],request.POST['num3_X2'],request.POST['num3_Y2'],request.POST['num3_RF']],
            [request.POST['num4_vol'],request.POST['num4_X1'],request.POST['num4_Y1'],request.POST['num4_X2'],request.POST['num4_Y2'],request.POST['num4_RF']],
            [request.POST['num5_vol'],request.POST['num5_X1'],request.POST['num5_Y1'],request.POST['num5_X2'],request.POST['num5_Y2'],request.POST['num5_RF']],
            [request.POST['num6_vol'],request.POST['num6_X1'],request.POST['num6_Y1'],request.POST['num6_X2'],request.POST['num6_Y2'],request.POST['num6_RF']],
            ],
            index=["1","2","3","4","5","6"],
            columns=["Vol","X1","Y1","X2","Y2","RF"])
        print(df1)
        df1.to_csv("Electrode-D5G.csv")
        #df1.to_csv("path/file_name.csv")
        #with open("Electrode-D5G.csv #csvの場所確認#", "w") as D5G:
        
        
       
    





        #ここにCSVに書き込む指令を出す
        return render(request, 'MD/ef.html', context)

class ecView(View):
    def get(self, request, *args, **kwargs):
        context = {
            'message': "実行してください",
        }
        return render(request, 'MD/ec.html', context)


    def post(self, request, *args, **kwargs):
        print("実行中...")
        context = {
            'message':"実行中...",
        }
        #ここにshell起動コマンド
        cmd = ["./field.sh", "0", "1", "100", "Erectrode-D5G.csv"]
        subprocess.check_call(cmd, shell=False)
        #subprocess.run(["field.sh","0", "1", "100", "Erectrode-D5G.csv"])#これだと開けない
        #os.system("field.sh 0 1 100 Erectrode-D5G.csv") #引数が渡せない ファイルが開いちゃう。
        print("しゅうりょう！おめでとう！！！！")
        context = {
            'message':"終了しました",
        }
        return render(request, 'MD/ec.html', context)













class MDcView(View):
    def get(self, request, *args, **kwargs):
        context = {
            'massage':"計算してください",
        }
        return render(request, 'MD/MDc.html', context)
    
    def post(self, request, *args, **kwargs):
        print("MDやるよ０－ん")
        context = {
            'massage':"実行中...",
        }
        #ここにpy起動コマンド
        #mpirun -n 8 python MDCylinderFinal.py D5G Particle-Stokes.csv Electrode-D5G.csv 1030 1000 50000 100 1000000 1 50 0
        command = [ "python.exe", "MDCylinderFinal.py", "D5G", "Particle-Stokes.csv", "Electrode-D5G.csv", "1030", "1000", "50000", "100", "1000000", "1", "50", "0", ]
        subprocess.run(command)


        print("otyukare <3")
        context = {
            'massage':"MD計算完了",
        }
        return render(request, 'MD/MDc.html', context)



ef = efView.as_view()
ec = ecView.as_view()
MDc = MDcView.as_view()