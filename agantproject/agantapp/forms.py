from django import forms

class SearchForm(forms.Form):
    pdbid = forms.CharField(label='Enter PDBID', max_length=10)