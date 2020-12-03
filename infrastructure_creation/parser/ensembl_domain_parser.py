"""
This script creates genome positions of each protein domains for each proteins.
This script parses the downloaded content by ensembl_proteinid_domain_fetch.py
"""


import os
import sys
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
# Check the directory for this script, create one if not exist.
data_repo_dir = create_dir(output_dir, OUTPUT_DATA_REPO)


def ensembl_domains(ensembl_domain_raw_data):
    """
    This is to find out which databases are already downloaded
    :param ensembl_domain_raw_data: Downloaded content by ensembl_proteinid_domain_fetch.py
    :return: Set of domains.
    """
    db = set()  # Initialize set to fill up
    for t in ensembl_domain_raw_data:  # For each protein in the downloaded content
        for i in ensembl_domain_raw_data[t]:  # For each domain entry of a given protein
            db.add(i["type"])  # Add the database name to the set
    return db


def gene3d_domain_parser(ensembl_domain_raw_data):
    """
    Fetches Gene3d data from a downloaded database.
    :param ensembl_domain_raw_data: Downloaded content by ensembl_proteinid_domain_fetch.py
    :return: Dictionary mapping to proteins to domain ids, domain ids to list of ranges.
    """
    output = dict()  # Initialize a dictionary to fill up with protein id keys

    for protein_id in ensembl_domain_raw_data:  # For each protein in the downloaded content
        gene3d_domains_ranges = dict()  # Create a dictionary to fill up with domain id keys

        for domain_entry in ensembl_domain_raw_data[protein_id]:  # For each domain in a given protein
            if domain_entry["type"] == 'Gene3D':  # If it is not a Gene3D entry, ignore.

                # Check if domain id is saved. Append the list with start and end positions (amino acid positions)
                if domain_entry["id"] not in gene3d_domains_ranges:
                    gene3d_domains_ranges[domain_entry["id"]] = [[domain_entry["start"], domain_entry["end"]]]
                else:
                    gene3d_domains_ranges[domain_entry["id"]].append([domain_entry["start"], domain_entry["end"]])

        if gene3d_domains_ranges:  # If gene3d_domains_ranges is populated at least one domain entry
            output[protein_id] = gene3d_domains_ranges  # Add the protein to the output dictionary
    # NOTES: This function cannot differentiate between:
    # - Two domains with the same domain id in a protein
    # - One domain composed of two segments in a protein
    # This function assumes they are the same.
    return output


def genome_range_to_positions(genome_ranges, reverse_strand):
    """
    Convert the ranges in protein_id_genome_map[protein_id]["gnPositions"] to genomic coordinates.
    This function is designed considering the format of the protein_id_genome_map[protein_id]["gnPositions"]
    :param genome_ranges: protein_id_genome_map[protein_id]["gnPositions"]. See protein_id_genome_map.py
    :param reverse_strand: Boolean to indicate whether the protein is on reverse orientation or forward
    :return: Numpy array of genomic coordinates. The order follows the order of amino acids.
    """
    # Note: i <= j always. But ranges in genome_ranges are in reverse orientation if reverse_strand = True.
    if reverse_strand:  # If reverse stand
        # [::-1] is to reverse the order of ranges before making the whole array as two dimensional.
        return np.concatenate([np.arange(i, j + 1) for i, j in genome_ranges][::-1])
    else:  # If forward strand
        return np.concatenate([np.arange(i, j + 1) for i, j in genome_ranges])


def domain_to_genome_positions(ensembl_domain_filtered, protein_id_genome_map):
    """
    As the name suggest, for each domain entry in ensembl_domain_filtered, calculate the corresponding genome positions.
    :param ensembl_domain_filtered: Dictionary, gene3d_domain_parser output
    :param protein_id_genome_map: Dictionary, see protein_id_genome_map.py
    :return: Dictionary for genomic ranges, chromosome information and orientation of the protein domains.
    """

    def helper(protein_ranges, gn_pos):
        """
        Determines the amino acid ranges' corresponding genomic positions.
        :param protein_ranges: Nested list. Each list contains amino acid start and end positions.
        :param gn_pos: np.array composed of genomic positions. gnPositions elements follows the amino acid order.
        :return:
        """
        genome_positions = [gn_pos[(j - 1) * 3: k * 3] for j, k in protein_ranges]  # *3 → 1 aa = 3 nt
        return np.concatenate(genome_positions)

    output = dict()  # Initialize a dictionary to fill up with protein id keys

    for ind, protein_id in enumerate(ensembl_domain_filtered):  # For each protein in the downloaded content
        progressBarForTerminal(ind, len(ensembl_domain_filtered) - 1)  # Print progress bar for current iteration

        chromosome = protein_id_genome_map[protein_id]["chromosome"]  # Save the chromosome name
        reverse_strand = protein_id_genome_map[protein_id]["reverse_strand"]  # Save the orientation
        # Create a entry in the output dictionary with the protein id as key.
        # Initialize gene3d dictionary in this dictionary
        output[protein_id] = {"chromosome": chromosome, "reverse_strand": reverse_strand, "gene3d": {}}

        # Get the genome positions for the protein from list of genome ranges. Takes into account the orientation.
        gn_positions = genome_range_to_positions(protein_id_genome_map[protein_id]["gnPositions"], reverse_strand)
        for domain_entry in ensembl_domain_filtered[protein_id]:  # For each entry domain entry
            # Calculate which genome positions the ranges correspond to. Add the entry under gene3d key.
            output[protein_id]["gene3d"][domain_entry] = helper(ensembl_domain_filtered[protein_id][domain_entry], gn_positions)

    return output


print(f"{bcolors.HEADER}Relevant Joblib objects are loading to RAM.{bcolors.ENDC}")
# Import the downloaded content. See ensembl_proteinid_domain_fetch.py
domain_map = joblib.load(os.path.join(output_dir, "protein_structure/proteinid_domains_map.joblib"))
genome_map = joblib.load(os.path.join(output_dir, "protein_structure/protein_id_genome_map.joblib"))

domains = ensembl_domains(domain_map)  # See which domains are there to report
print(f"Available domains in downloaded content: {', '.join(domains)}")  # Print available domains
print("Only Gene3D is implemented to the analysis pipeline.")

# Parse Gene3D domains
print(f"{bcolors.HEADER}Gene3D domains are being parsed.{bcolors.ENDC}")
domain_map = gene3d_domain_parser(domain_map)  # Parse the downloaded content
# Keep only proteins which have entries in both genome map and domain map.
genome_map, domain_map = sync_dictionaries(genome_map, domain_map)

print(f"{bcolors.HEADER}Domain ranges are converted into genomic coordinates.{bcolors.ENDC}")
# Calculate genome coordinates for the protein ranges of domains
ensembl_gene3d_genome_coordinates = domain_to_genome_positions(domain_map, genome_map)


# Write down the output dictionary and list as Joblib object for convenience in later uses.
print(f"{bcolors.HEADER}Results are written to disk: {data_repo_dir}{bcolors.ENDC}")
os.chdir(data_repo_dir)  # Since all outputs will be written in the same directory.
joblib.dump(ensembl_gene3d_genome_coordinates, "ensembl_gene3d_genome_coordinates.joblib")
print("Done: ensembl_gene3d_genome_coordinates.joblib")


# End of the script
