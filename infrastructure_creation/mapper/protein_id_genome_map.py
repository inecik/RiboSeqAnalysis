"""
This script is to create a dictionary mapping amino acid positions in ensembl protein to genomic coordinates.
"""

import os
import sys
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
# Check the directory for this script, create one if not exist.
data_repo_dir = create_dir(output_dir, OUTPUT_DATA_REPO)
temp_repo_dir = create_dir(temp_dir, TEMP_DATA_REPO)


print(f"{bcolors.HEADER}Transcript of genes are being fetched from MANE project.{bcolors.ENDC}")
ero = ensembl_release_object_creator()  # Create a new ensembl release object of pyensembl library.
# Below lines are the same as the script for downloading.
# The aim is to get the same set of protein_ids here as well.
# See ensembl_proteinid_domain_fetch.py script for explanation of below lines.
gene_map_mane = mane_gene_maps_with_defined_protein(ero, temp_repo_dir)
transcript_uniprot_map = db_to_uniprot_mapping(temp_repo_dir, mapping="Ensembl_TRS")
uniprot_protein_map = uniprot_to_db_mapping_multiple(temp_repo_dir, "Ensembl_PRO")
uniprot_list = list(set([k for i in gene_map_mane
                         for j in ero.transcript_ids_of_gene_id(i)
                         if transcript_uniprot_map.get(j)
                         for k in transcript_uniprot_map[j]]))
protein_ids = list(set([j for i in uniprot_list for j in uniprot_protein_map[i]["Ensembl_PRO"]]))
print(f"{bcolors.HEADER}Unnecessary temporary variables are being deleted from RAM.{bcolors.ENDC}")
del gene_map_mane, transcript_uniprot_map, uniprot_protein_map, uniprot_list


# Main calculation of protein to genomic coordinates starts from here
print(f"{bcolors.HEADER}Protein ranges are being converted to genome coordinates.{bcolors.ENDC}")
counter1 = 0  # Initialize a counter
counter2 = 0  # Initialize a counter
protein_id_genome_map = dict()  # Initialize a dictionary to fill up

for ind, protein_id in enumerate(protein_ids):  # Iterate through protein IDs
    progressBarForTerminal(ind, len(protein_ids) - 1)  # Print progress bar to report the progress
    try:
        # Find associated transcript if exist. Raises ValueError if no transcript is found.
        t = ero.transcript_by_protein_id(protein_id)
    except ValueError:
        continue  # Ignore protein if no transcripts of it found.

    if not t.is_protein_coding:  # Ignore protein if its transcript is not labelled as protein coding
        continue

    # Get the CDS ranges of the transcript, convert list of tuples to list of lists.
    cds_ranges = [list(i) for i in t.coding_sequence_position_ranges]  # Note: Required for changing end position below

    protein_length = len(t.protein_sequence)  # Calculate the length of protein sequence
    cds_length = range_sum(cds_ranges)  # Calculate the total length of CDS sequences
    assert t.strand in ["-", "+"]  # Make sure transcript's strand is either "-" or "+"
    reverse_strand = t.strand == "-"  # Create a 'reverse_strand' variable, True if transcript is on reverse strand.

    if protein_length * 3 == cds_length:  # If protein length and CDS length seems consistent
        protein_id_genome_map[protein_id] = {  # Add the CDS ranges to dictionary
            "gnPositions": cds_ranges,
            "chromosome": t.contig,
            "reverse_strand": reverse_strand
        }
        continue  # Ignore the remaining part of the loop. Move to the next protein

    # Below lines are created for some of the proteins, which has inconsistency between protein length and CDS length.
    if t.contains_start_codon and t.contains_stop_codon:
        counter1 += 1  # Ignore if both start and stop exist but there is a inconsistency anyway.


    elif t.contains_start_codon and not t.contains_stop_codon:  # If there is a start codon but no stop codon
        # Below line was different initially. When I realized no transcript is under that condition, I deleted.
        # No example of it seen (theoretically possible though).
        assert protein_length * 3 <= cds_length  # Added here to make sure
        # Extend CDS to allow protein_length * 3 == cds_length
        # Removing nucleotides from end of CDS will not shift the frame for mapping protein to genome coordinates
        if not reverse_strand:  # If the transcript is on positive strand
            if cds_ranges[-1][1] - (cds_length - protein_length * 3) >= cds_ranges[-1][0]:
                # Above is to make sure after removal, start < end
                cds_ranges[-1][1] -= cds_length - protein_length * 3
            else:
                raise Exception(f"{bcolors.FAIL}Error: Unexpectedly short exon{bcolors.ENDC}")
        else: # If the transcript is on negative strand
            if cds_ranges[-1][0] - (cds_length - protein_length * 3) <= cds_ranges[-1][1]:
                # Above is to make sure after removal, start < end
                cds_ranges[-1][0] += cds_length - protein_length * 3
            else:
                raise Exception(f"{bcolors.FAIL}Error: Unexpectedly short exon{bcolors.ENDC}")

    # The calculation of it is not easy, can create frame shift mistakes.
    elif not t.contains_start_codon:
        counter2 += 1  # Ignore if no start codon exist

    # Check if protein length and CDS length seems consistent after removal of nucleotides from CDS
    if protein_length * 3 == range_sum(cds_ranges):
        protein_id_genome_map[protein_id] = {  # Add the CDS ranges to dictionary
            "gnPositions": cds_ranges,
            "chromosome": t.contig,
            "reverse_strand": reverse_strand
        }

# Print out some the ignored cases.
print("CDS length does not match with protein length:")
print(f"Contains both start and stop codon: {counter1}")
print(f"Does not contain start codon: {counter2}")


# Write down the output dictionary and list as Joblib object for convenience in later uses.
print(f"{bcolors.HEADER}Results are written to disk: {data_repo_dir}{bcolors.ENDC}")
os.chdir(data_repo_dir)  # Since all outputs will be written in the same directory.
joblib.dump(protein_id_genome_map, "protein_id_genome_map.joblib")
print("Done: protein_id_genome_map.joblib")


# End of script
