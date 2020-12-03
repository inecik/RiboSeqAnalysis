"""
This script creates genome positions of each protein domains for each proteins.
This script parses the downloaded content by uniprot_protein_structures_and_genome_positions_fetch.py
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


def filter_uniprot_features(raw_feature_script_output):
    """
    Filters raw output of the downloaded content. Passes relevant information only.
    Filtering is not perfect. It removes a lot of useful information as well.
    Important data such as nuclear localization signal is lost under "Motif".
    Only data which can be easily quantified is kept.
    :param raw_feature_script_output: One output file in uniprot_protein_structures_and_genome_positions_fetch.py
    :return:
    """

    relevant_info = dict()  # Initialize a dictionary to fill up with protein id keys

    for protein_id in raw_feature_script_output:  # For each entry in downloaded content. i is protein_id
        features = list()  # Initialize a list to fill up with protein id keys

        for feature_entry in raw_feature_script_output[protein_id]["features"]:  # For each feature entry for a protein
            if feature_entry["category"] not in (  # Filter if the entry is in following 'categories'.
                    "MUTAGENESIS",
                    "VARIANTS",
            ):
                if feature_entry["type"] not in (  # Filter if the entry is in following 'types'.
                        "NON_TER", "CONFLICT", "NON_STD",  # under category: SEQUENCE_INFORMATION
                        "DOMAIN", "REGION", "REPEAT", "BINDING", "SITE", "MOTIF",  # under category: DOMAINS_AND_SITES
                        "CHAIN", "PEPTIDE",  # under category: MOLECULE_PROCESSING
                ):
                    features.append(feature_entry)  # If non of above applies, add to the list

        if features:  # If list is populated at least once
            relevant_info[protein_id] = {"features": features,  # add the list in the dictionary
                                         "sequence": raw_feature_script_output[protein_id]["sequence"]}
            # Also keep the sequence information to use in a later section in the pipeline
    return relevant_info


def filter_correct_length_transcript(filtered_uniprot_to_genome_map, filtered_uniprot_features):
    """
    Filters out transcripts with CDS which does not consistent with corresponding protein's length
    :param filtered_uniprot_to_genome_map: Dictionary. The format is similar with output of uniprot_to_genome_map.py
    :param filtered_uniprot_features: Dictionary. The format is similar with feature_extract_to_range function output
    :return: Dictionary. Filtered version of the filtered_uniprot_to_genome_map
    """
    counter = 0  # Initialize a counter to report the number of filtered entry at the end
    output = dict()  # Initialize a dictionary to fill up with uniprot protein id as keys

    for p in filtered_uniprot_to_genome_map:  # for uniprot id in the genome_map
        # Get amino acid sequence and calculate CDS length by multiplying with three
        other_dict_value = len(filtered_uniprot_features[p]["sequence"]) * 3
        # Get transcripts of the defined under uniprot protein
        transcripts = list(filtered_uniprot_to_genome_map[p].keys())
        transcript_infos = list()  # Create a empty list to use as filtering

        for t in transcripts:  # for each transcript id defined under uniprot protein
            # Check if CDS of it is at expected length
            if other_dict_value != len(filtered_uniprot_to_genome_map[p][t]["gnPositions"]):
                transcript_infos.append(None)  # If not add None to the transcript_infos list
                counter += 1
            else:
                transcript_infos.append(filtered_uniprot_to_genome_map[p][t])  # Otherwise, append the content

        if any(transcript_infos):  # if transcript_infos is populated with at least one entry
            transcripts = list(compress(transcripts, transcript_infos))  # Remove wrong length transcript
            transcript_infos = list(compress(transcript_infos, transcript_infos))  # Remove wrong length transcript info
            # Add the Uniprot ID under output dict, with a dictionary of valid transcripts
            output[p] = {transcripts[i]: transcript_infos[i] for i in range(len(transcripts))}

    if counter:  # If at least one transcript deleted, report the number of transcript deleted
        print(f"Number of transcripts deleted: {counter}")

    return output


def define_uniprot_annotation_tracks(uniprot_features_dict):
    """
    This script is to define which annotation tracks (features) is going to be included to the pipeline.
    :param uniprot_features_dict: Dictionary, filter_uniprot_features output
    :return: Two sets composed of strings.
    """

    tracks = set()  # Initialize a set to fill up with annotation tracks
    tracks_compbias = set()  # Initialize a set to fill up with annotation tracks (of a specific type; compbias entries)

    for protein_id in uniprot_features_dict:  # For each protein in the downloaded and filtered content
        for feature_entry in uniprot_features_dict[protein_id]["features"]:  # For each feature in a single protein

            # If a feature among one of the following type
            if feature_entry["type"] in [
                "ZN_FING", "NP_BIND", "DNA_BIND", "DISULFID",
                "PROPEP", "TRANSMEM", "INTRAMEM", "ACT_SITE",
                "METAL", "MOD_RES", "LIPID", "CARBOHYD",
                "CROSSLNK", "HELIX", "STRAND","TURN",
                "TOPO_DOM", "CA_BIND", "SIGNAL", "TRANSIT", "INIT_MET"
            ]:
                tracks.add(feature_entry["type"])  # Add the feature entry under tracks set

            else:  # Otherwise
                if feature_entry["type"] == "COMPBIAS":  # If the entry is a compbias type entry
                    # Remove everything after '-rich' and add to the tracks_compbias set.
                    # Removal of '-rich' allows to end up with less tracks_compbias length, looks cleaner.
                    tracks_compbias.add(feature_entry["description"].split("-rich")[0])
                else:
                    # If 'description' is not a defined
                    if "description" not in feature_entry:
                        # Add the feature entry under tracks set in addition to the above list.
                        tracks.add(feature_entry["type"])  # It is actually just to make sure
                    else:  # Raise an error otherwise
                        raise Exception(f"{bcolors.FAIL}Error in define_uniprot_annotation_tracks.\n"
                                        f"A type with description is not evaluated.{bcolors.ENDC}")
    return tracks, tracks_compbias


def feature_extract_to_range(protein_features, tracks, compbias):
    """
    Creates nested lists for tracks and compbias.
    Nested lists shows where start and end positions of protein features.
    :param protein_features: Dictionary, filter_uniprot_features output
    :param tracks: Protein annotations. One of the output of define_uniprot_annotation_tracks
    :param compbias: Protein annotations. One of the output of define_uniprot_annotation_tracks
    :return: Updated version of protein_features with ranges of annotations
    """
    for protein_id in protein_features:  # for each protein id in the dictionary

        annotation = dict()  # Initialize a dict to fill up with annotation tracks
        annotation_compbias = dict()  # Initialize a dict to fill up with annotation tracks

        protein_info = protein_features[protein_id]  # Get the protein information from the dictionary
        for domain_entry in protein_info["features"]:  # For a domain entry of the protein

            if domain_entry["type"] == "COMPBIAS":  # If the feature type is "COMPBIAS"
                entry_compbias = domain_entry["description"].split("-rich")[0]  # Remove after '-rich' as before
                assert entry_compbias in compbias  # Make sure the feature is a defined compbias type

                try:
                    begin = int(domain_entry["begin"].replace("~", ""))  # Remove approximate sign
                    end = int(domain_entry["end"].replace("~", ""))  # Remove approximate sign
                except ValueError:  # ValueError if int function cannot supplied with numbers
                    # Note: Sometimes only approximate sign exists. It occurs when it is not known exactly.
                    continue  # Ignore those

                # Update annotation_compbias entry with start and end positions.
                if entry_compbias not in annotation_compbias:
                    annotation_compbias[entry_compbias] = [[begin, end]]
                else:
                    annotation_compbias[entry_compbias].append([begin, end])

            elif domain_entry["type"] in tracks:  # If the feature type is in 'tracks'

                try:  # Similar to above, remove approximate sign.
                    begin = int(domain_entry["begin"].replace("~", ""))
                except ValueError:  # ValueError if int function cannot supplied with numbers
                    continue

                try:  # Similar to above, remove approximate sign.
                    end = int(domain_entry["end"].replace("~", ""))
                except ValueError:  # ValueError if int function cannot supplied with numbers
                    # Only for "TRANSIT", "SIGNAL". I assumed that end position is 25 more of begin position
                    if domain_entry["type"] in ["TRANSIT", "SIGNAL"] and begin < 5:  # If begin position is less than 5
                        end = begin + 25  # 25 arbitrarily. I considered it is an average value after looking the data
                    else:
                        continue  # Ignore otherwise

                # Update annotation entry with start and end positions.
                if domain_entry["type"] not in annotation:
                    annotation[domain_entry["type"]] = [[begin, end]]
                else:
                    annotation[domain_entry["type"]].append([begin, end])

            else:  # Unexpected condition. Raise error.
                raise Exception(f"{bcolors.FAIL}Error in feature_extract. Unexpected feature type found.{bcolors.ENDC}")

        del protein_features[protein_id]["features"]  # Remove the features dictionary
        protein_features[protein_id]["tracks"] = annotation  # Add tracks as the key and annotation as the value
        protein_features[protein_id]["compbias"] = annotation_compbias  # Similar to above line, update the dictionary

    return protein_features


def feature_ranges_to_genome(protein_features_ranges, filtered_uniprot_to_genome_map):
    """
    As the name suggest, for each domain entry in protein_features_ranges, calculate the corresponding genome positions.
    :param protein_features_ranges: Dictionary. The format is similar with feature_extract_to_range output
    :param filtered_uniprot_to_genome_map: Dictionary. The format is similar with output of uniprot_to_genome_map.py
    :return: Dictionary for genomic ranges, chromosome information and orientation of the protein domains.
    """

    def helper(gn_positions, protein_ranges):
        """
        Determines the amino acid ranges' corresponding genomic positions.
        :param gn_positions: Genomic positions of the transcript, which is in the same order with amino acid sequence
        :param protein_ranges: Nested list. List composed of start and end amino acid positions in protein
        :return: numpy array of genomic positions of the range
        """
        genome_positions = list()  # Initialize a list to populate
        for j, k in protein_ranges:
            genome_positions.append(gn_positions[(j - 1) * 3: k * 3])  # Get the ranges' corresponding genome positions
        return np.concatenate(genome_positions)  # Concatenate nested list to have a two dimensional np array

    output = dict()  # Initialize a dictionary to fill up with protein id keys

    for ind, protein_id in enumerate(protein_features_ranges):  # For each protein in the downloaded content
        progressBarForTerminal(ind, len(protein_features_ranges) - 1)  # Print progress bar for current iteration

        tracks = protein_features_ranges[protein_id]["tracks"]  # Get the 'tracks' features
        compbias = protein_features_ranges[protein_id]["compbias"]  # Get the 'compbias' features
        # Initialize protein_id in output dictionary. Initialize transcripts in output[protein_id]
        output[protein_id] = {i: None for i in filtered_uniprot_to_genome_map[protein_id]}

        for transcript_id in filtered_uniprot_to_genome_map[protein_id]:  # for each transcript associated with protein
            # Get the genomic positions of the transcript, which is in the same order with amino acid sequence
            position_info = filtered_uniprot_to_genome_map[protein_id][transcript_id]["gnPositions"]
            # For each feature ranges in tracks and compbias, get the genomic positions
            tracks_positions = {i: helper(position_info, tracks[i]) for i in tracks}
            compbias_positions = {i: helper(position_info, compbias[i]) for i in compbias}
            # Add resulting dictionaries together with chromosome and orientation
            output[protein_id][transcript_id] = {
                "chromosome": filtered_uniprot_to_genome_map[protein_id][transcript_id]["chromosome"],
                "reverse_strand": filtered_uniprot_to_genome_map[protein_id][transcript_id]["reverse_strand"],
                "tracks": tracks_positions,
                "compbias": compbias_positions,
            }
    return output


print(f"{bcolors.HEADER}Relevant Joblib objects are loading to RAM.{bcolors.ENDC}")
# Import the downloaded content. See uniprot_protein_structures_and_genome_positions_fetch.py
feature_map = joblib.load(os.path.join(output_dir, "protein_structure/uniprotid_features_map.joblib"))
# Filter out unnecessary entries in order to reduce the memory consumption and simply later steps.
feature_map = filter_uniprot_features(feature_map)
# Determine which tracks and compbias (as annotation tracks) will be used.
ann_tracks, ann_compbias = define_uniprot_annotation_tracks(feature_map)
# Get the selected features' protein ranges for the dowloaded and filtered content
feature_map = feature_extract_to_range(feature_map, ann_tracks, ann_compbias)
# Import the created variable. See uniprot_to_genome_map.py
genome_map = joblib.load(os.path.join(output_dir, "protein_structure/uniprot_to_genome_map.joblib"))
# Keep only proteins which have entries in both genome map and domain map.
feature_map, genome_map = sync_dictionaries(feature_map, genome_map)
# Remove transcripts which are not at expected length by comparing the amino acid length with CDS length of transcripts
genome_map = filter_correct_length_transcript(genome_map, feature_map)
# Keep only proteins which have entries in both genome map and domain map after the filtering
feature_map, genome_map = sync_dictionaries(feature_map, genome_map)


print(f"{bcolors.HEADER}Feature ranges are converted into genomic coordinates.{bcolors.ENDC}")
uniprot_feature_genome_coordinates = feature_ranges_to_genome(feature_map, genome_map) # todo


# Write down the output dictionary and list as Joblib object for convenience in later uses.
print(f"{bcolors.HEADER}Results are written to disk: {data_repo_dir}{bcolors.ENDC}")
os.chdir(data_repo_dir)  # Since all outputs will be written in the same directory.
joblib.dump(uniprot_feature_genome_coordinates, "uniprot_feature_genome_coordinates.joblib")
print("Done: uniprot_feature_genome_coordinates.joblib")


# End of the script
