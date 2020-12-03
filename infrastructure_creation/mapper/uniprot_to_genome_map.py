"""
This script is to create a dictionary mapping amino acid positions in Uniprot protein to genomic coordinates.
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


print(f"{bcolors.HEADER}Relevant Joblib objects are loading to RAM.{bcolors.ENDC}")
# Import the downloaded content. See uniprot_protein_structures_and_genome_positions_fetch.py
uniprot_to_genome_map = joblib.load(os.path.join(output_dir, "protein_structure/uniprotid_positions_map.joblib"))


def genomic_location_helper(entry):
    """
    This is to get a numpy array of genomic positions. The order of numpy array follows the amino acid positions.
    :param entry: JSON response of the server. It is for a single Uniprot protein.
    :return: dictionary of chromosome (str), reverse_strand (bool), gnPositions (np.array)
    """
    chromosome = entry["genomicLocation"]["chromosome"]  # Save the chromosome of Uniprot entry
    reverse_strand = entry["genomicLocation"]["reverseStrand"]  # Save the strand of Uniprot entry
    genomic_ranges = list()  # Initialize a list to fill up in the following loop

    # Below is to determine CDS ranges of a protein
    for exon in entry["genomicLocation"]["exon"]:  # For each exon entry in the response
        info = exon["genomeLocation"]  # Save the genomeLocation entry of the exon

        if len(info) == 1 and "position" in info:  # If there is only one in info, no start and stop position defined.
            g_start = g_end = info["position"]["position"]  # Set both start and stop the same

        elif len(info) == 2:  # In case two in info
            # Make sure start and stop positions are there
            assert "begin" in info and "end" in info, f"No end or begin in entry:\n{info}"
            g_start = info["begin"]["position"]  # Save start position
            g_end = info["end"]["position"]  # Save stop position
        else:  # No example of it is seen but added here to make sure
            raise Exception(f"{bcolors.FAIL}Unexpected error #1:\n{info}{bcolors.ENDC}")
        genomic_ranges.append((g_start, g_end))  # Add to the list

    # Below is to create a numpy array from CDS ranges which is in the same direction with protein.
    plus_or_minus = -1 if reverse_strand else 1  # Determine the strand. -1 if reverse, +1 if forward orientation
    # Calculate the genomic positions by taking into account the orientation.
    gnPositions = np.concatenate([np.arange(i, j + plus_or_minus, plus_or_minus) for i, j in genomic_ranges])
    # Do not sort to make sure 'exon' (actually CDS) order is represented correctly in gnPositions.
    mapping = {"chromosome": chromosome,
               "reverse_strand": reverse_strand,
               "gnPositions": gnPositions}
    # Create a dictionary and report the resulting genomic positions.
    return mapping


# Main calculation of protein to genomic coordinates starts from here
print(f"{bcolors.HEADER}Protein ranges are being converted to genome coordinates.{bcolors.ENDC}")
for ind, p in enumerate(uniprot_to_genome_map):  # For each uniprot response from the server.
    progressBarForTerminal(ind, len(uniprot_to_genome_map) - 1)  # Print progress bar to report the progress

    for i in uniprot_to_genome_map[p]["gnCoordinate"]:  # For each transcript entry for a given Uniprot protein
        assert i['ensemblTranscriptId'] not in uniprot_to_genome_map[p]  # Check if the same transcript entries exist
        # Parse the genomic positions by genomic_location_helper
        # Add the results to the dictionary with the key corresponding transcript-id
        # Do not create a dictionary use existing one (Updates the server response actually)
        uniprot_to_genome_map[p][i['ensemblTranscriptId']] = genomic_location_helper(i)

    # Delete every other information not required from the response. (Removes unnecessary information from response)
    for i in ["gnCoordinate", "gene", "name", "taxid", "protein", "accession", "sequence"]:
        if i in uniprot_to_genome_map[p]:
            del uniprot_to_genome_map[p][i]


# Write down the output dictionary and list as Joblib object for convenience in later uses.
print(f"{bcolors.HEADER}Results are written to disk: {data_repo_dir}{bcolors.ENDC}")
os.chdir(data_repo_dir)  # Since all outputs will be written in the same directory.
joblib.dump(uniprot_to_genome_map, "uniprot_to_genome_map.joblib")
print("Done: uniprot_to_genome_map.joblib")


# End of script
