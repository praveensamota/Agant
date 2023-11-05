import csv
from django.core.management.base import BaseCommand
from agantapp.models import ProteinData

class Command(BaseCommand):
    help = 'Import protein data from CSV files'

    def handle(self, *args, **options):
        # Open and read antagonist CSV file
        with open('/home/praveen/agant/agantproject/agantapp/management/commands/pdb_titles_ant.csv', 'r') as antagonist_file:
            antagonist_reader = csv.DictReader(antagonist_file)
            for row in antagonist_reader:
                self.save_protein_data(row)

        # Open and read agonist CSV file
        with open('/home/praveen/agant/agantproject/agantapp/management/commands/pdb_titles_ag.csv', 'r') as agonist_file:
            agonist_reader = csv.DictReader(agonist_file)
            for row in agonist_reader:
                self.save_protein_data(row)

        self.stdout.write(self.style.SUCCESS('Data imported successfully.'))
        
    def save_protein_data(self, row):
        # Handle empty values and convert to appropriate data types
        lig_vol = float(row['Lig_vol']) if row['Lig_vol'] else 0.0
        pock_vol = float(row['Pock_vol']) if row['Pock_vol'] else 0.0
        drug_score = float(row['DrugScore']) if row['DrugScore'] else 0.0
        average_b_factor = float(row['Average B-Factor']) if row['Average B-Factor'] else 0.0
        all_atoms_rsa_cavity = float(row['All_atoms_rsa_cavity']) if row['All_atoms_rsa_cavity'] else 0.0
        predicted_average_pkd = float(row['Predicted Average pKd']) if row['Predicted Average pKd'] else 0.0
        total_surface_area = float(row['Total Surface Area (A^2)']) if row['Total Surface Area (A^2)'] else 0.0

        protein_data = ProteinData(
            PDB=row['PDB'],
            Title=row['Title'],
            Response=row['Response'],
            Ligand=row['Ligand'],
            Uniprot_accession=row['Uniprot_accession'],
            Lig_vol=lig_vol,
            Pock_vol=pock_vol,
            DrugScore=drug_score,
            Average_B_Factor=average_b_factor,
            All_atoms_rsa_cavity=all_atoms_rsa_cavity,
            Predicted_Average_pKd=predicted_average_pkd,
            Total_Surface_Area=total_surface_area
        )
        protein_data.save()

