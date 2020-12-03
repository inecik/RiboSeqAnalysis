"""
This script is to fetch data from the Uniprot database using for structural information of proteins.
It uses already existed API but it does not parse the resulting JSON. 
The script is designed to not crash when the Internet connection is lost.
It saves the resulting dictionary as Joblib object.
"""


import sys
import os
import requests
import joblib
from random import randint
from time import sleep, strftime

repository_name = "RiboSeqAnalysis"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from infrastructure_creation.functions import *
from archieve.common_functions import *


# Inputs
output_dir = sys.argv[1]
temp_dir = sys.argv[2]


# File Names
OUTPUT_DATA_REPO = "protein_structure"
TEMP_DATA_REPO = "protein_structure"
# Check the directory for this script, create one if not exist.
data_repo_dir = create_dir(output_dir, OUTPUT_DATA_REPO)
temp_repo_dir = create_dir(temp_dir, TEMP_DATA_REPO)


print(f"{bcolors.HEADER}Transcript of genes are being fetched from MANE project.{bcolors.ENDC}")
ero = ensembl_release_object_creator()  # Create a new ensembl release object of pyensembl library.
# Determine which genes we are interested in so that download only relevant proteins.
gene_map_mane = mane_gene_maps_with_defined_protein(ero, temp_repo_dir)


print(f"{bcolors.HEADER}Helper variables are being creating.{bcolors.ENDC}") # To determine which Uniprot IDs are of interest
# Create a dictionary which maps ensembl transcript id to uniprot id
transcript_uniprot_map = db_to_uniprot_mapping(temp_repo_dir, mapping="Ensembl_TRS")
# Create a unique list of uniprot IDs for which we download the relevant structural information
uniprot_list = list(set([k for i in gene_map_mane  # for each gene in the MANE project
                         for j in ero.transcript_ids_of_gene_id(i)  # Find out the transcripts of the genes
                         if transcript_uniprot_map.get(j)  # If transcript has associated Uniprot ID.
                         for k in transcript_uniprot_map[j]]))  # Get the Uniprot ID of the transcript
# Whole procedure ensures all proteins are from a gene in MANE project, and has an associated Uniprot ID.


print(f"{bcolors.HEADER}Unnecessary temporary variables are being deleted from RAM.{bcolors.ENDC}")
# We are not going to use these variables in the remaining part of the script anyway. 
del gene_map_mane, ero, transcript_uniprot_map


print(f"{bcolors.HEADER}Structural informations are being fetched from EBI REST server.{bcolors.ENDC}")
# Following two is to fetch data about proteins structural information
uniprotid_features_map = dict()  # Initialize dictionary to save downloaded content
uniprotid_features_map_unfound = list()  # Initialize list to save proteins for which no information could be downloaded.
# Following two is to fetch the data to help mapping proteins on genome 
uniprotid_positions_map = dict()  # Initialize dictionary to save downloaded content
uniprotid_positions_map_unfound = list()  # Initialize list to save proteins for which no information could be downloaded.


iteration = 0  # Initialize variable to use in the following while loop
disconnected = False  # Initialize variable to use in the following while loop

while iteration != len(uniprot_list):  # Continue until having response from the server for all proteins.
    uniprot_id = uniprot_list[iteration]  # Protein of interest in current iteration 
    
    try:
        # Raises ConnectionError, TimeoutError, or ReadTimeout if a problem occurs in the Internet connection
        query2 = requests.get(f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}",
                              headers={"Accept": "application/json"}, timeout=5)
        query3 = requests.get(f"https://www.ebi.ac.uk/proteins/api/coordinates/{uniprot_id}",
                              headers={"Accept": "application/json"}, timeout=5)
        # Raise TimeoutError if no response is received. Continue if responses are received for both query. 
        
        # For feature information
        if query2.status_code != 200:  # If a problem found in the response
            uniprotid_features_map_unfound.append((uniprot_id, query2.status_code))  # Update the relevant list
        else:  # If no problem found
            uniprotid_features_map[uniprot_id] = query2.json()  # Add the response to the relevant dictionary
        
        # For genome positions
        if query3.status_code != 200:  # If a problem found in the response
            uniprotid_positions_map_unfound.append((uniprot_id, query3.status_code))  # Update the relevant list
        else:  # If no problem found
            uniprotid_positions_map[uniprot_id] = query3.json()  # Add the response to the relevant dictionary
        
        if disconnected:  # Print following if when returning from a connection related error.
            print(f"Reconnected:     {strftime('%Y-%m-%d %H:%M:%S')}")
            disconnected = False
        
        # If code can reach here without any error, do following..
        progressBarForTerminal(iteration, len(uniprot_list) - 1,  # Print progress bar for current iteration
                               prefix =f"{bcolors.OKCYAN}Downloading:", barLength=30,
                               suffix =f"{bcolors.BOLD}    {iteration+1}/{len(uniprot_list)}{bcolors.ENDC}")
        iteration += 1  # To move on the next element in the list of proteins

    except (requests.exceptions.ConnectionError, TimeoutError, requests.exceptions.ReadTimeout):
        
        if not disconnected:  # Print following if at the first time during connection.
            print(f"\nConnection lost: {strftime('%Y-%m-%d %H:%M:%S')}")
            disconnected = True
        
        # Try to request information from the server again after 2-7 seconds later.
        sleep(randint(2, 7))


# Write down the output dictionary and list as Joblib object for convenience in later uses.
print(f"{bcolors.HEADER}Results are written to disk: {data_repo_dir}{bcolors.ENDC}")
os.chdir(data_repo_dir)  # Since all outputs will be written in the same directory.

joblib.dump(uniprotid_positions_map_unfound, "uniprotid_positions_map_unfound.joblib")
print("Done: uniprotid_positions_map_unfound.joblib")

joblib.dump(uniprotid_positions_map, "uniprotid_positions_map.joblib")
print("Done: uniprotid_positions_map.joblib")

joblib.dump(uniprotid_features_map_unfound, "uniprotid_features_map_unfound.joblib")
print("Done: uniprotid_features_map_unfound.joblib")

joblib.dump(uniprotid_features_map, "uniprotid_features_map.joblib")
print("Done: uniprotid_features_map.joblib")


# End of the script
