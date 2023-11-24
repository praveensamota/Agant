from django.shortcuts import render
from .forms import SearchForm
from .models import agant
from agantapp.models import ProteinData
from django.core.paginator import Paginator
from django.http import JsonResponse
from django.db.models import Q
from django.urls import reverse
from .forms import PredictionForm
from django.core.files.storage import FileSystemStorage
from agantapp.utils.ml_scripts import rootfunction
import os
import pandas as pd
import pickle
import subprocess



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

current_directory = 'agantapp'



# # Load the trained XGBoost model from the .pkl file
# model_filename_pkl = os.path.join(current_directory, 'xgboost_model.pkl')
# with open(model_filename_pkl, 'rb') as file:
#     loaded_model_pkl = pickle.load(file)
# print("Trained XGBoost model loaded successfully from .pkl file")

# Assuming you have a Django form for file uploads (PredictionForm)
def prediction_view(request):
    if request.method == 'POST':
        form = PredictionForm(request.POST, request.FILES)
        if form.is_valid():
            # Get user inputs from the form
            pdb_file = request.FILES['pdb_file']
            ligand_code = form.cleaned_data['ligand_code']
            SMILES = form.cleaned_data['SMILES']
            mol2_file = request.FILES['mol2_file']
            pdb_id = form.cleaned_data['pdb_id']

            # Save the uploaded files to a temporary location

            temp_directory = os.path.join(current_directory, 'temp')
            print(temp_directory)
            if not os.path.exists(temp_directory):
                os.makedirs(temp_directory)
            pdb_file_path = os.path.join(current_directory, 'temp', pdb_file.name)
            print("PDB File Path:", pdb_file_path)
            with open(pdb_file_path, 'wb') as destination:
                for chunk in pdb_file.chunks():
                    destination.write(chunk)

            mol2_file_path = os.path.join(current_directory, 'temp', mol2_file.name)
            with open(mol2_file_path, 'wb') as destination:
                for chunk in mol2_file.chunks():
                    destination.write(chunk)

            # Run the necessary scripts to process the input files
            output_directory = os.path.join(current_directory, 'temp_output')
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)

            # Print file paths for debugging
            
            print("MOL2 File Path:", mol2_file_path)

            predictions=rootfunction(pdb_file,ligand_code,SMILES,mol2_file,pdb_id,current_directory,output_directory)

            # # Run dpocket
            # run_dpocket(pdb_id, ligand_code, output_directory)

            # # Generate PaDel descriptors
            # generate_descriptors(SMILES, output_directory)

            # # Run NACCESS on the cavity
            # selected_cavity_file = os.path.join(output_directory, f'{pdb_id}_cavity_1.pdb')
            # naccess_output_file_cavity, csv_output_file_cavity = process_pdb_file(selected_cavity_file)

            # # Calculate percentages of helices and sheets
            # ss_percentages = calculate_percentages(selected_cavity_file)

            # # Process the selected surface file for H-bond interactions
            # selected_surface_file = os.path.join(output_directory, f'{pdb_id}_surface_1.pdb')
            # output_csv_hbond = os.path.join(output_directory, f'H-BOND_interactions_surface.csv')
            # process_single_pdb(selected_surface_file, ligand_code, output_csv_hbond)

            # # Load the features for prediction
            # Trained_features = os.path.join(current_directory, 'trained_features.txt')
            # with open(Trained_features, 'r') as file:
            #     feature_list = [line.strip() for line in file]

            # # Combine CSV files
            # output_combined_csv = os.path.join(output_directory, 'combined_data_all_files.csv')
            # os.system(f"python combine_csv.py {output_directory} {output_combined_csv}")

            # # Load the combined CSV file
            # combined_csv = pd.read_csv(output_combined_csv)

            # # Select only the columns corresponding to the features in the list
            # selected_data = combined_csv[feature_list]

            # # Make predictions using the loaded XGBoost model
            # predictions_pkl = loaded_model_pkl.predict(selected_data)

            # # Display or use the predictions as needed
            # print("Predictions (from .pkl):", predictions_pkl)

            # Prepare the context to be passed to the template
            context = {
                'predictions': predictions,
                # 'ss_percentages': ss_percentages,
                # 'output_csv_hbond': output_csv_hbond,
            }

            return render(request, 'agantapp/prediction_template.html', context)

    else:
        form = PredictionForm()

    return render(request, 'agantapp/prediction_template.html', {'form': form})
