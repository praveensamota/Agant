from django.urls import path
from . import views 

urlpatterns = [
    path("", views.index_view, name="index"),
    path('antagonist/', views.antagonist_view, name='antagonist_view'),
    path('agonist/', views.agonist_view, name='agonist_view'),
    path('search/', views.live_search_view, name='live_search'),
    path('prediction/', views.prediction_view, name='prediction'),
]