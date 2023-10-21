from django.shortcuts import render
from .forms import SearchForm
from .models import agant
# Create your views here.

def index_view(request):
    if request.method == 'POST':
        form = SearchForm(request.POST)
        if form.is_valid():
            pdbid = form.cleaned_data['pdbid']
            agent = agant.objects.filter(pdbid=pdbid).first()
            if agent:
                context = {
                    'form': form,
                    'agent': agent,
                    'agent_type': agent._type,
                }
            else:
                context = {
                    'form': form,
                    'error_message': 'No matching record found for the given PDBID.',
                }
        else:
            context = {
                'form': form,
            }
    else:
        form = SearchForm()
        context = {
            'form': form,
        }
    return render(request, 'agantapp/index.html', context)
