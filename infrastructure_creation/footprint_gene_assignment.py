"""
This script contains functions to calculate the number of footprints on genes.
The output for a gene is an array of integers with the same length of the gene length.
It saves the resulting dictionary as Joblib object.
"""

import os
import sys
import re
import pysam
import numpy as np
import joblib

repository_name = "RiboSeqAnalysis"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from infrastructure_creation.functions import *
from archieve.common_functions import *


def footprint_assignment(sam_path, assignment=-15, verbose=False):
    """
    This is to assign the footprints into chromosomes and then into genomic positions.
    :param verbose: Boolean. Report or not report the progress.
    :param sam_path: UMI deduplicated SAM file to use for genome assignment
    :param assignment: Integer to count from an end. Negative numbers are assignment from 3'. '-15' → Arpat et. al 2020
    :return: Dictionary of assigned positions
    """
    # If 'verbose' activated, count the number of lines in SAM file to be able to print a progress bar.
    if verbose:
        with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open SAM file with pysam library
            iteration = sam_handle.count()  # Call count method to find out total read count.
    
    with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open SAM file with pysam library
        sam_iterator = sam_handle.fetch()  # Get the iterator to iterate through reads in for loop
        assigned_positions = dict()  # Initialize dictionary with keys composed of chromosome names
        
        for ind, e in enumerate(sam_iterator):  # Iterate through the entries
            
            # If 'verbose' activated, print the progress bar in a meaningful intervals
            if verbose and (ind % 1000 == 0 or ind == iteration - 1):
                progressBarForTerminal(ind, iteration - 1, barLength=20)

            assert sum(e.get_cigar_stats()[1][5:]) == 0  # Make sure cigar string composed of D, I, M, N, S only
            # other than 'D', get_reference_positions method perfectly works.
            # See the following link to understand why we should keep 'D'
            # https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/cigar.md

            if 'D' not in e.cigarstring:
                reference_positions = e.get_reference_positions()  # Gets the genomic coordinates of the aligned region
            else:
                # If there is a matching element (M) before to 'D', merge D to M.
                for regex_match in re.finditer(r"([0-9]*)M[0-9]*D", e.cigarstring):
                    # regex_match.group() is like '48M2D'
                    m_len, d_len = re.findall(r"[0-9]+", regex_match.group())  # Get the numbers
                    m_len = int(m_len) + int(d_len)
                    e.cigarstring = e.cigarstring.replace(regex_match.group(), f"{m_len}M")
                assert not re.search(r"([0-9]*D)", e.cigarstring)  # Make sure D's are merged to adjacent M
                reference_positions = e.get_reference_positions()  # Get genomic coordinates including 'D'

            if e.reference_name not in assigned_positions:  # If the chromosome name is not present in the dictionary
                assigned_positions[e.reference_name] = [reference_positions[assignment]]
            else:  # Add the assigned position of the footprint into the relevant list in the dictionary
                assigned_positions[e.reference_name].append(reference_positions[assignment])

    for chromosome in assigned_positions:  # When assignment ends for each entry in the chromosome dictionary
        assigned_positions[chromosome] = np.sort(np.array(assigned_positions[chromosome]))  # Make it np.array and sort

    return assigned_positions


def chromosomes_genes_matcher(ensembl_release_object, gene_list):
    """
    This script is to assign genes into chromosomes.
    :param ensembl_release_object: Created by ensembl_release_object_creator() function
    :param gene_list: List of genes to assign chromosomes
    :return: Dictionary of assigned genes.
    """
    chromosomes = dict()  # Initialize a dictionary to fill up
    for gene in gene_list:  # For each gene in the gene_list
        geob = ensembl_release_object.gene_by_id(gene)  # Create an instance of Ensembl Gene object
        # Find the chromosome name (contig), add the gene into the relevant list
        chromosomes[geob.contig] = chromosomes.get(geob.contig, []) + [geob.id]
    return chromosomes


def gene_to_array(ensembl_release_object, gene_id):
    """
    This script is to create a array of genomic positions for a gene.
    :param ensembl_release_object: Created by ensembl_release_object_creator() function
    :param gene_id: String, ensembl gene id
    :return: Numpy array of genomic coordinates
    """
    geob = ensembl_release_object.gene_by_id(gene_id)  # Create an instance of Ensembl Gene object
    assert geob.start < geob.end  # Make sure start position is smaller than end position.
    # Note: Even genes in reverse strand, geob.start is smaller than geob.end
    return np.arange(geob.start, geob.end + 1)  # +1 → both ends should be included


def footprint_counts_to_genes(footprint_genome_assignment, chromosome_gene_map, positions_gene_map, verbose=False):
    """
    Assigns the footprints to genes. Output shows how many footprints are assigned at each position of a gene.
    :param footprint_genome_assignment: Output of footprint_assignment function
    :param chromosome_gene_map: Output of chromosomes_genes_matcher function
    :param positions_gene_map: Dictionary mapping gene_id to array of genomic coordinates, output of gene_to_array func.
    :param verbose: Boolean. Report or not report the progress.
    :return: Dictionary mapping gene_id to numpy array (with the same length of gene), which compose of integers.
    """
    gene_footprint_assignment = dict()  # Initialize a dictionary to fill up

    for ind, chr_name in enumerate(footprint_genome_assignment):  # For each chromosome
        if verbose:  # If 'verbose' activated, print the progress bar
            progressBarForTerminal(ind, len(footprint_genome_assignment) - 1,
                                   suffix=f"    Chromosome: {chr_name}", barLength=20)
        # Get all genomic coordinates assigned to a footprint
        chr_pos_footprints = footprint_genome_assignment[chr_name]
        # Get the number of unique elements and the counts of these unique elements
        element_unique, element_count = np.unique(chr_pos_footprints, return_counts=True)
        # Create a dictionary from this information, mapping genomic position with the number of footprints assigned
        footprints_counts = dict(zip(element_unique, element_count))

        for gene_id in chromosome_gene_map[chr_name]:  # For each gene in the same chromosome
            pos_genes = positions_gene_map[gene_id]  # Get the genomic positions of the gene
            # For each position in the gene, look at the footprints_counts, get the number if found, 0 otherwise.
            gene_footprint_assignment[gene_id] = np.array([footprints_counts.get(i, 0) for i in pos_genes])
    return gene_footprint_assignment


def gene_assignment_main(identifier, sam_path, output_directory, temp_directory, verbose=True):
    """
    Controls the gene assignment procedure.
    :param identifier: Name of the SAM file in order to append to the joblib object's name
    :param sam_path: UMI deduplicated SAM file to use for genome assignment.
    :param output_directory: Output directory, see main.py
    :param temp_directory: Temp directory, see main.py
    :param verbose: Boolean. Report or not report the progress.
    :return: None. It creates a Joblib object.
    """

    # File Names
    output_data_repo = "protein_structure"
    temp_data_repo = "protein_structure"
    # Check temp directory for this script, create one if not exist.
    data_repo_dir = create_dir(output_directory, output_data_repo)
    temp_repo_dir = create_dir(temp_directory, temp_data_repo)

    if verbose:
        print(f"{bcolors.HEADER}Genes are being fetched from MANE project.{bcolors.ENDC}")
    ero = ensembl_release_object_creator()  # Create a new ensembl release object of pyensembl library.
    # Determine which genes we are interested in so that assign footprints to these genes only.
    gene_list = list(mane_gene_maps_with_defined_protein(ero, temp_repo_dir).keys())
    # Get the mapping from chromosomes to list of genes which they carry
    chromosome_gene = chromosomes_genes_matcher(ero, gene_list)
    # Get array of genomic positions for all genes, create a dictionary out of it.
    positions_gene = {i: gene_to_array(ero, i) for i in gene_list}

    if verbose:
        print(f"{bcolors.HEADER}Footprint are being assigned to genomic coordinates.{bcolors.ENDC}")
    # Assign footprints to genomic positions
    footprint_genome_assignment = footprint_assignment(sam_path, verbose=verbose)
    # Remove chromosomes which does not contain footprint assignment, or which does not contain genes
    chromosome_gene, footprint_genome_assignment = sync_dictionaries(chromosome_gene, footprint_genome_assignment)

    if verbose:
        print(f"{bcolors.HEADER}Footprint counts are being calculated and assigned to genes.{bcolors.ENDC}")
    # Assign footprints to gene positions
    footprints_on_genes = footprint_counts_to_genes(footprint_genome_assignment, chromosome_gene, positions_gene, verbose=verbose)

    # Write down the output dictionary and list as Joblib object for convenience in later uses.
    if verbose:
        print(f"{bcolors.HEADER}Results are written to disk: {data_repo_dir}{bcolors.ENDC}")
    os.chdir(data_repo_dir)  # Since all outputs will be written in the same directory.
    joblib.dump(footprints_on_genes, f"footprints_on_genes_{identifier}.joblib")
    if verbose:
        print(f"Done: footprints_on_genes_{identifier}.joblib")


if __name__ == '__main__':
    assert len(sys.argv) == 5  # Make sure correct number of arguments are supplied
    sam_name = sys.argv[1]  # Append this string to the output path
    sam_file = sys.argv[2]  # Absolute path of sam file
    output_dir = sys.argv[3]  # Absolute path of output directory
    temp_dir = sys.argv[4]  # Absolute path of temp directory
    gene_assignment_main(sam_name, sam_file, output_dir, temp_dir)  # Call the main function


# End of the script
