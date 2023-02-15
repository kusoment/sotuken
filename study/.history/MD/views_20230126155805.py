from django.shortcuts import render
from django.http import HttpResponse

# Create your views here.
def index(request):
    context = {
        'foo': 'foo string!',
    }
    return render(request, 'MD/index.html', context)

def ef(request):
    return render(request, "MD/ef.html",)