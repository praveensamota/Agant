from django.shortcuts import render
from .forms import SearchForm
from .models import agant
from agantapp.models import ProteinData
from django.core.paginator import Paginator
from django.http import JsonResponse
from django.db.models import Q
from django.urls import reverse
# Create your views here.

def index_view(request):
    return render(request, 'agantapp/index.html')


def antagonist_view(request):
    queryset = ProteinData.objects.filter(Response='Antagonist')  # Replace this with your queryset
    paginator = Paginator(queryset, 150)  # Paginate by 150 items per page

    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)

    return render(request, 'agantapp/antagonist_template.html', {'page_obj': page_obj})

def live_search_view(request):
    query = request.GET.get('query', '')
    response_type = request.GET.get('type', '')
    results = ProteinData.objects.filter(
        Q(Response=response_type) & (Q(PDB__icontains=query) | Q(Title__icontains=query))
    )
    data = []
    for result in results:
        entry_data = {
            'PDB': result.PDB,
            'Title': result.Title,
            'Response': result.Response,
            'Ligand': result.Ligand,
            'Uniprot_accession': result.Uniprot_accession,
            'Lig_vol': result.Lig_vol,
            'Pock_vol': result.Pock_vol,
            'DrugScore': result.DrugScore,
            'Average_B_Factor': result.Average_B_Factor,
            'All_atoms_rsa_cavity': result.All_atoms_rsa_cavity,
            'Predicted_Average_pKd': result.Predicted_Average_pKd,
            'Total_Surface_Area': result.Total_Surface_Area,
            'PDB_URL': f'https://www.rcsb.org/structure/{result.PDB}',
        }
        data.append(entry_data)
    return JsonResponse(data, safe=False)

def agonist_view(request):
    agonist_data = ProteinData.objects.filter(Response='Agonist')  # Replace 'Agonist' with your actual value for Response field
    paginator = Paginator(agonist_data, 150)  # Paginate by 150 items per page

    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)

    return render(request, 'agantapp/agonist_template.html', {'page_obj': page_obj})