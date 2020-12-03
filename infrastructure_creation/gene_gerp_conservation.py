"""
This script fetches GERP scores for genes.
"""

import os
import sys
import pyBigWig
import numpy as np
import joblib

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
GERP_BIGWIG = "gerp_conservation_scores.homo_sapiens.GRCh38.bw"
# Check the directory for this script, create one if not exist.
data_repo_dir = create_dir(output_dir, OUTPUT_DATA_REPO)
temp_repo_dir = create_dir(temp_dir, TEMP_DATA_REPO)


# Download conservation score if not exists already
gerp_bigwig_path = os.path.join(temp_repo_dir, GERP_BIGWIG)  # Path for gerp conservation score containing bigwig file
if not os.access(gerp_bigwig_path, os.R_OK) or not os.path.isfile(gerp_bigwig_path):  # If not exist already.
    print(f"{bcolors.HEADER}Human GERP scores is being downloading from the server.{bcolors.ENDC}")
    subprocess.run((f"cd {temp_repo_dir}; "  # Change directory to save the file
                   f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/compara/conservation_scores/"
                   f"88_mammals.gerp_conservation_score/{GERP_BIGWIG}"), shell=True)  # Curl the database


print(f"{bcolors.HEADER}Genes are being fetched from MANE project.{bcolors.ENDC}")
ero = ensembl_release_object_creator()  # Create a new ensembl release object of pyensembl library.
# Determine which genes we are interested in so that fetch data for only relevant ones.
gene_list = list(mane_gene_maps_with_defined_protein(ero, temp_repo_dir).keys())


print(f"{bcolors.HEADER}Conservation scores are being fetched from BigWig file.{bcolors.ENDC}")
gene_gerp = dict()  # Initialize a dictionary to fill up with protein id keys

with pyBigWig.open(gerp_bigwig_path) as bw:  # Open BigWig file with pyBigWig library
    gerp_chroms = list(bw.chroms().keys())  # Get the names of chromosomes in the BigWig file

    for ind, g in enumerate(gene_list):  # Iterate over the genes
        progressBarForTerminal(ind, len(gene_list) - 1)  # Print the progress bar

        gene = ero.gene_by_id(g)  # Create a gene object to get the start and end positions of it
        if gene.contig not in gerp_chroms:  # If gene's chromosome not in BigWig file
            # Create an array of np.nan values with the same length of gene
            gene_gerp[gene.id] = np.full(gene.length, np.nan)
        else:
            # Create an array of conservation score values, which is at the same length of gene
            gene_gerp[gene.id] = bw.values(gene.contig, gene.start, gene.end + 1, numpy=True)


# Write down the output dictionary and list as Joblib object for convenience in later uses.
print(f"{bcolors.HEADER}Results are written to disk: {data_repo_dir}{bcolors.ENDC}")
os.chdir(data_repo_dir)  # Since all outputs will be written in the same directory.
joblib.dump(gene_gerp, "genes_positions_cons-scores.joblib")
print("Done: genes_positions_cons-scores.joblib")


# End of the script
