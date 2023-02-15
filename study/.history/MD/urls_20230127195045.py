from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('ef/', views.ef, name="ef"),
    path('ec/', views.ec, name="ec"),
    path('MDc/', views.MDc, name="MDc"),
    path('plot/', views.plot, name="plot"),
    path('set_plot', views.get_svg, name="set_plot")
]