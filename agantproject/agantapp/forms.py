from django import forms

class SearchForm(forms.Form):
    pdbid = forms.CharField(label='Enter PDBID', max_length=10)


class PredictionForm(forms.Form):
    pdb_file = forms.FileField()
    pdb_id = forms.CharField()
    ligand_code = forms.CharField()
    SMILES = forms.CharField(max_length=255, required=True)
    mol2_file = forms.FileField()