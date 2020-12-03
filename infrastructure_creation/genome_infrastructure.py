"""
This script creates a dictionary for each gene by combining all sorts of relevant information.
Structural information, evolutionary conservation information, transcript ranges, CDS ranges, footprint assignment.
"""

import os
import sys
import numpy as np
import joblib
from functools import reduce

repository_name = "RiboSeqAnalysis"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from infrastructure_creation.functions import *
from archieve.common_functions import *



# Temp and output directories
output_dir = sys.argv[1]
temp_dir = sys.argv[2]


# File Names and directory creation
OUTPUT_DATA_REPO = "protein_structure"
TEMP_DATA_REPO = "protein_structure"
# Check the directory for this script, create one if not exist.
data_repo_dir = create_dir(output_dir, OUTPUT_DATA_REPO)
temp_repo_dir = create_dir(temp_dir, TEMP_DATA_REPO)


# Load previously created variables
print(f"{bcolors.HEADER}Relevant Joblib objects are loading to RAM.{bcolors.ENDC}")
domains_gene3d = joblib.load(os.path.join(data_repo_dir, "ensembl_gene3d_genome_coordinates.joblib"))
features_uniprot = joblib.load(os.path.join(data_repo_dir, "uniprot_feature_genome_coordinates.joblib"))
gene_conservation = joblib.load(os.path.join(data_repo_dir, "genes_positions_cons-scores.joblib"))
footprints_on_genes = joblib.load(os.path.join(data_repo_dir, "footprints_on_genes.joblib"))


# Create some dictionaries to use in the loop
print(f"{bcolors.HEADER}Helper variables are being creating.{bcolors.ENDC}")
ero = ensembl_release_object_creator()  # Create a new ensembl release object of pyensembl library.
# Create a mapping from transcript id to ensembl protein ID
enst_to_ensp = {ero.transcript_id_of_protein_id(i): i for i in ero.protein_ids()}
# Genes are being fetched from MANE project
gene_map_mane = mane_gene_maps_with_defined_protein(ero, temp_repo_dir)
# Create a mapping from transcript id to uniprot protein id
enst_to_uniprot = {}  # Initialize a dictionary
for k, v in features_uniprot.items():  # for each protein id and associated transcript names in features_uniprot
    for i in v:  # for each transcript in transcripts
        enst_to_uniprot[i] = enst_to_uniprot.get(i, []) + [k]  # Add the mapping to the dictionary


def range_genome_mapper(genome_positions, genome_ranges):
    """
    A helper function to create a boolean array to show which positions are covered by genome_ranges
    :param genome_positions: np array of genome coordinates
    :param genome_ranges: Nested list of ranges
    :return: Boolean numpy array
    """
    for rgm_j_test, rgm_k_test in genome_ranges:  # Make sure first element is smaller or equal to the second
        assert rgm_j_test <= rgm_k_test
    range_positions = np.concatenate([np.arange(rgm_j, rgm_k + 1) for rgm_j, rgm_k in genome_ranges])
    return np.isin(genome_positions, range_positions)


counter_chromosome_strand_unmatch = 0  # Initialize a counter to report inconsistencies
genome_base = dict()  # Initialize a dictionary to fill up with gene information
feature_index = dict()  # Initialize a dictionary to map features with genes in which they are found
gene_ids = list(gene_map_mane.keys())  # Get the list of gene names

for ind, gene_id in enumerate(gene_ids):  # For each gene in the list
    progressBarForTerminal(ind, len(gene_ids) - 1, suffix=f"    {ind+1} / {len(gene_ids)}", barLength=25)  # Progress

    geOb = ero.gene_by_id(gene_id)  # Create a gene object to get some information about the gene
    assert geOb.start < geOb.end  # Make sure start is always smaller than end position, regarless of orientation
    # Get transcripts of the genes which are protein coding (have a protein ensembl id)
    valid_transcripts = [i for i in ero.transcript_ids_of_gene_id(geOb.id) if i in enst_to_ensp]
    valid_proteins = [enst_to_ensp[i] for i in valid_transcripts]  # Get associated protein id
    positions = np.arange(geOb.start, geOb.end + 1)  # Calculate genomic coordinates of gene

    gene_page = {  # Create a new dictionary to fill up some information about the gene
        "gene_id": geOb.id,  # Gene ID
        "gene_name": geOb.gene_name,  # Gene name
        "contig": geOb.contig,  # Chromosome
        "strand": geOb.strand,  # '+' or '-' strand
        "start": geOb.start,  # Start position (including)
        "end": geOb.end,  # End position (including)
        "transcript_ids": valid_transcripts,  # Protein coding transcripts
        "protein_ids": valid_proteins,  # Protein coding transcripts' proteins
        "mane_transcript": gene_map_mane[geOb.id],  # Transcript ID chosen by MANE project
        # Get boolean array of CDS positions for MANE transcript
        "mane_transcript_cds_PA": range_genome_mapper(positions, ero.transcript_by_id(
            gene_map_mane[geOb.id]).coding_sequence_position_ranges),
        # Get boolean array of exon positions for all exons the gene contains
        "exon_ranges_PA": range_genome_mapper(positions, gene_exon_ranges(ero, geOb.id)),
        # Get boolean array of CDS positions for all transcripts CDSs of the gene
        "cds_ranges_PA": range_genome_mapper(positions, gene_cds_ranges(ero, geOb.id)),
        # Find the transcripts with longest CDS, return boolean array of CDS positions
        "longest_cds_with_defined_protein_PA": range_genome_mapper(positions, gene_longest_cds_with_defined_protein(
            ero, geOb.id)),
        "conservation": gene_conservation[geOb.id],  # Get the evolutionary conservation scores of the gene positions
        "footprints": footprints_on_genes[geOb.id]  # Get the footprint assignment scores for the gene positions
    }

    # Gene3d annotations
    gene3d = dict()  # Initialize a dictionary to fill up with gene3d information
    for valid_protein in valid_proteins:  # For each protein of the gene
        if valid_protein in domains_gene3d:  # If it has a annotation by Ensembl for Gene3D domains

            try:
                # Test whether chromosome and orientation is consistent
                assert domains_gene3d[valid_protein]["chromosome"] == gene_page["contig"]
                assert domains_gene3d[valid_protein]["reverse_strand"] == (gene_page["strand"] == "-")
            except AssertionError:
                counter_chromosome_strand_unmatch += 1  # If not consistent report ignored cases.
                continue  # If not consistent ignore the proein

            for gene3d_id in domains_gene3d[valid_protein]["gene3d"]:  # for each defined domain for the protein
                p_dom_info = domains_gene3d[valid_protein]["gene3d"][gene3d_id]  # get the information
                gene3d[gene3d_id] = gene3d.get(gene3d_id, []) + [p_dom_info]  # add to the dictionary with its id as key

    for gene3d_reduce in gene3d:  # Get the union of all annotations for a given domain
        # It is required since many proteins of a gene exists so domain annotation is redundant for a region.
        # It ensures all annotation for a gene are obtained somehow. However, I should test the method for its validity.
        reduced_gene3d = reduce(np.union1d, gene3d[gene3d_reduce])  # Get the union of a domain entries
        gene3d[gene3d_reduce] = np.isin(positions, reduced_gene3d)  # Convert the genome positions to boolean array

    # Uniprot annotations
    compbias = dict()  # Initialize a dictionary to fill up with compbias information
    tracks = dict()  # Initialize a dictionary to fill up with tracks information
    for valid_transcript in valid_transcripts:  # For each transcript of the gene
        for uniprot_id in enst_to_uniprot.get(valid_transcript, []):  # For Uniprot ID of the transcript if exists

            try:
                # Test whether chromosome and orientation is consistent
                assert features_uniprot[uniprot_id][valid_transcript]["chromosome"] == gene_page["contig"]
                assert features_uniprot[uniprot_id][valid_transcript]["reverse_strand"] == (gene_page["strand"] == "-")
            except AssertionError:
                counter_chromosome_strand_unmatch += 1  # If not consistent report ignored cases.
                continue  # If not consistent ignore the proein

            for t_compbias in features_uniprot[uniprot_id][valid_transcript]["compbias"]:  # For each defined compbias
                # Get the information and add to the dictionary with its id as the key
                t_comp_info = features_uniprot[uniprot_id][valid_transcript]["compbias"][t_compbias]
                compbias[t_compbias] = compbias.get(t_compbias, []) + [t_comp_info]

            for t_tracks in features_uniprot[uniprot_id][valid_transcript]["tracks"]:  # For each defined tracks
                # Get the information and add to the dictionary with its id as the key
                t_track_info = features_uniprot[uniprot_id][valid_transcript]["tracks"][t_tracks]
                tracks[t_tracks] = tracks.get(t_tracks, []) + [t_track_info]

    for compbias_reduce in compbias:  # Get the union of all annotations for a given feature
        reduced_compbias = reduce(np.union1d, compbias[compbias_reduce])  # Get the union of a domain entries
        compbias[compbias_reduce] = np.isin(positions, reduced_compbias)  # Convert the genome positions to boolean arr.

    for tracks_reduce in tracks:  # Get the union of all annotations for a given feature
        reduced_tracks = reduce(np.union1d, tracks[tracks_reduce])  # Get the union of a domain entries
        tracks[tracks_reduce] = np.isin(positions, reduced_tracks)  # Convert the genome positions to boolean array

    # Add the infos to the gene dictionary
    gene_page["gene3d"] = gene3d
    gene_page["compbias"] = compbias
    gene_page["tracks"] = tracks

    # Get the features or domains for this gene. Save them to the dictionary for later analysis
    for i in gene3d.keys() | compbias.keys() | tracks.keys():
        feature_index[i] = feature_index.get(i, []) + [geOb.id]  # Feature/domain to list of genes mapping

    # Save to the main dictionary
    genome_base[geOb.id] = gene_page  # Save the gene dictionary to the output dictionary, which consists of all genes.

if counter_chromosome_strand_unmatch > 0:  # If there is a consistency, print following as warning.
    print(f"{bcolors.WARNING}counter_chromosome_strand_unmatch: {counter_chromosome_strand_unmatch}{bcolors.ENDC}")


# Write down the output dictionary and list as Joblib object for convenience in later uses.
print(f"{bcolors.HEADER}Results are written to disk JOBLIB: {data_repo_dir}{bcolors.ENDC}")
os.chdir(data_repo_dir)

joblib.dump(genome_base, "genome_base.joblib")
print("Done: genome_base.joblib")

joblib.dump(feature_index, "genome_base_feature_index.joblib")
print("Done: genome_base_feature_index.joblib")


# End of the script
