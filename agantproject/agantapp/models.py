from django.db import models

# Create your models here.
class agant(models.Model):
    pdbid = models.CharField(max_length=200)
    l_id = models.CharField(max_length=100)
    _type = models.CharField(max_length=100)
    uniprotid = models.URLField(max_length=100)
    
    def __str__(self):
        return self.pdbid

class ProteinData(models.Model):
    PDB = models.CharField(max_length=100)
    Title = models.CharField(max_length=255)
    Response = models.CharField(max_length=255)
    Ligand = models.CharField(max_length=255)
    Uniprot_accession = models.CharField(max_length=255)
    Lig_vol = models.FloatField(null=True)
    Pock_vol = models.FloatField(null=True)
    DrugScore = models.FloatField(null=True)
    Average_B_Factor = models.FloatField(null=True)
    All_atoms_rsa_cavity = models.FloatField(null=True)
    Predicted_Average_pKd = models.FloatField(null=True)
    Total_Surface_Area = models.FloatField(null=True)

    def __str__(self):
        return self.PDB
