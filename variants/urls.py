from . import views
from django.urls import path

urlpatterns = [
    path('', views.index, name='index'),
    path('predict/<str:login>', views.predict, name='predict'),
    path('mutationMap/<str:rsid>/<str:login>', views.mutationMap, name='mutationMap'), #<str:rsid>
    path('help', views.help, name='help'),
    path('about', views.about, name='about'),
]