from django.shortcuts import render
from django.http import HttpResponse

# Create your views here.
def index(request):
    context = {
        'foo': 'foo string!',
    }
    return render(reuest, 'MD/index.html', context)