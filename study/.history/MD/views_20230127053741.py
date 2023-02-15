
from django.shortcuts import render
from django.http import HttpResponse
from django.views import View
import csv

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
        #print(request.POST['num1_vol'],request.POST['num2_vol'])
        with open("Electrode-D5G.csv #csvの場所確認#", "w") as D5G:
            open
        
       
    





        #ここにCSVに書き込む指令を出す
            return render(request, 'MD/ef.html', context)

ef = efView.as_view()