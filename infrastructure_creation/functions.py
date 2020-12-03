"""
This script contains supplemental functions to be used in this module.
"""

import gzip
import os
import subprocess
import sys
from itertools import compress
from shutil import which
from pyensembl import EnsemblRelease

repository_name = "RiboSeqAnalysis"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from archieve.common_functions import *


def uniprot_mapping_database_download(temp_repo_dir):
    """
    Function is to download or find the file in temp_repo_dir.
    :param temp_repo_dir: String. Directory to download or find the file
    :return: String. Path of dat file
    """
    # Download the data from the server
    mapping_info = "HUMAN_9606_idmapping.dat.gz"  # File name in the FTP server
    # Paths to download the database or look for the database if it already exists
    data_path_compressed = os.path.join(temp_repo_dir, mapping_info)  # Compressed .gz format
    data_path = os.path.splitext(data_path_compressed)[0]  # Uncompressed format
    if not os.access(data_path, os.R_OK) or not os.path.isfile(data_path):  # Check if the file exist
        subprocess.run((f"cd {temp_repo_dir}; curl -L -R -O ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/"
                        f"knowledgebase/idmapping/by_organism/{mapping_info}"), shell=True)  # Download otherwise
        subprocess.run(f"cd {temp_repo_dir}; gzip -d {data_path_compressed}", shell=True)  # Uncompress the file
    return data_path


def uniprot_to_db_mapping_multiple(temp_repo_dir, *args):
    """
    Parse Uniprot ID mapping to create a mapping between protein id to database of interest
    :param temp_repo_dir: String. Directory to download or find the file
    :param args: Strings such as "Ensembl", "Ensembl_TRS", "Ensembl_PRO", "Gene_Name", "Gene_Synonym"
    :return: Dictionary of Uniprot ID to dictionary ("Ensembl", "Ensembl_TRS") which consist of list of string items.
    """
    data_path = uniprot_mapping_database_download(temp_repo_dir)  # Download data or get the absolute path of it.
    with open(data_path, "r") as input_handle:  # Open the file
        data_raw = input_handle.readlines()  # Read its content
    df = dict()  # Create a dictionary to fill up in the following for loop
    for entry in data_raw:  # for each entry in the data HUMAN_9606_idmapping.dat
        entry = entry.strip().split("\t")  # Clean the line from flanking white spaces, split the line into cells.
        entry[0] = entry[0].split("-")[0]  # Remove '-' in Uniprot entries if exists
        if entry[1] in args:  # Check if the entry is among the requested ones
            # Add the protein entry appropriately.
            if entry[0] not in df:
                df[entry[0]] = {entry[1]: [entry[2]]}  # Add if uniprot ID is not exist in the output dictionary
            elif entry[1] not in df[entry[0]]:
                df[entry[0]][entry[1]] = [entry[2]]  # Add the database type under the uniprot ID
            else:
                df[entry[0]][entry[1]].append(entry[2])  # Add the entry only type under the database type
    return df


def db_to_uniprot_mapping(temp_repo_dir, mapping):
    """
    This is to get mapping dictionary to a database of interest to Uniprot (for Human)
    :param temp_repo_dir: String. Absolute path of directory to save necessary files
    :param mapping: String. Database of interest. Check the downloaded file to see the options
    :return: Dictionary of mapping.
    """
    data_path = uniprot_mapping_database_download(temp_repo_dir)  # Get the database file path
    with open(data_path, "r") as input_handle:  # Open the file
        data_raw = input_handle.readlines()  # Read its content
    df = dict()  # Initialize the dictionary
    for entry in data_raw:  # For each entry in the database
        entry = entry.strip().split("\t")  # Remove flanking white space, split the line into cells
        if entry[1] == mapping:  # If the line contains the information of interest
            entry[0] = entry[0].split("-")[0]  # Remove dash and afterwards in Uniprot ID
            if entry[2] not in df:  # Check if the entry is added to the mapping before
                df[entry[2]] = [entry[0]]
            else:
                df[entry[2]].append(entry[0])
    return df


def ensembl_release_object_creator(version=100):
    """
    This function is to create a ensembl release object or abort the script if the environment is not set up properly.
    :param version: Which version to use. The last version is 100.
    :return: Ensembl release object
    """
    try:
        return EnsemblRelease(version)
    except ValueError:
        sys.exit(f"Run following command before the code:\npyensembl install --release {version} --species human\n")


def genes_proteins(ensembl_release_object):
    """
    This function to create gene → protein map.
    :param ensembl_release_object: Created by ensembl_release_object_creator() function
    :return: Dictionary for gene → protein mapping
    """
    proteins = ensembl_release_object.protein_ids()  # Get all the proteins in the database
    gene_protein_map = dict()  # Create an empty dictionary to fill up
    for protein_id in proteins:  # For each protein
        gene_object = ensembl_release_object.gene_by_protein_id(protein_id)  # Get the gene id corresponding to protein
        if gene_object.id not in gene_protein_map:  # Look if the gene is already in the list
            gene_protein_map[gene_object.id] = [protein_id]  # Add list of protein_id as value, gene id as key
        else:  # Otherwise
            gene_protein_map[gene_object.id].append(protein_id)  # Append the list of protein ids
    return gene_protein_map


def range_sum(ranges):
    """
    The function is to calculate the total CDS length of a transcript.
    :param ranges: List of ranges. Obtained by coding_sequence_position_ranges method of transcript object
    :return: Integer. Total length of the CDS
    """
    total = 0  # Initialize
    for i, j in ranges:  # For each CDS of a transcript, get start and end positions
        # j always larger or equal to i
        total += (j - i + 1)  # Plus 1 is due to the fact that start and end is included.
    return total


def longest_cds_gene(ensembl_relase_object):
    """
    It reports the transcript with longest CDS disregarding biological function, being protein coding or not.
    :param ensembl_relase_object: Created by ensembl_relase_object_creator() function
    :return: Dictionary for gene → transcript mapping
    """
    gene_map = genes_proteins(ensembl_relase_object)  # Get gene to protein mapping
    for gene_id in gene_map:  # For each gene in the dictionary
        lengths = list()  # Initialize a list for CDS lengths
        transcripts = list()  # Initialize a list for transcript IDs
        for protein_id in gene_map[gene_id]:  # For each protein for a given gene
            transcript_object = ensembl_relase_object.transcript_by_protein_id(protein_id)  # Get the transcript ID
            lengths.append(range_sum(transcript_object.coding_sequence_position_ranges))  # Calculate the length of it
            transcripts.append(transcript_object.transcript_id)  # Append to the list
        gene_map[gene_id] = transcripts[lengths.index(max(lengths))]  # Report the longest one among the transcripts
    return gene_map


def longest_cds_gene_with_defined_protein(ensembl_relase_object, temp_repo_dir):
    """
    This function is similar to longest_cds_gene but this one takes into account if transcript has a Uniprot ID.
    :param ensembl_relase_object: Created by ensembl_relase_object_creator() function
    :param temp_repo_dir: String. Absolute path of directory to save necessary files
    :return: Dictionary for gene → transcript mapping
    """
    transcript_uniprot = db_to_uniprot_mapping(temp_repo_dir, mapping="Ensembl_TRS")  # Get transcript to Uniprot map
    gene_protein_map = genes_proteins(ensembl_relase_object)  # Get gene to protein mapping
    gene_transcript_map = dict()  # Initialize dictionary to report at the end
    for gene_id in gene_protein_map:  # For each gene in the dictionary
        lengths = list()  # Initialize a list for CDS lengths
        transcripts = list()  # Initialize a list for transcript IDs
        for protein_id in gene_protein_map[gene_id]:  # For each protein for a given gene
            transcript_object = ensembl_relase_object.transcript_by_protein_id(protein_id)  # Get the transcript ID
            lengths.append(range_sum(transcript_object.coding_sequence_position_ranges))  # Calculate the length of it
            transcripts.append(transcript_object.transcript_id)  # Append to the list
        # Create bool list to check if transcripts have a corresponding Uniprot entry
        compress_with = [t in transcript_uniprot for t in transcripts]
        transcripts = list(compress(transcripts, compress_with))  # Filter transcripts list based on the boolean filter
        lengths = list(compress(lengths, compress_with))  # Filter lengths list based on the boolean filter
        try:
            max_index = lengths.index(max(lengths))  # Get the index with maximum length
        except ValueError:  # If there is not transcript left after the filtering.
            continue  # Ignore the gene as it is obviously not a protein coding and no Uniprot ID
        gene_transcript_map[gene_id] = transcripts[max_index]  # Add transcript id to the mapping dictionary
    return gene_transcript_map


def mane_gene_transcript_map(temp_repo_dir):
    """
    Downloads and parses MANE initiative's data. See the link: https://www.ncbi.nlm.nih.gov/refseq/MANE/
    It only selects one transcript for a gene.
    MANE Select is chosen only, Plus Clinic is just ignored.
    :param temp_repo_dir: String. Absolute path of directory to save necessary files
    :return: Dictionary of Gene → Transcript.
    """
    mapping_info = "MANE.GRCh38.v0.92.summary.txt.gz"  # File name in the FTP server
    # Paths to download the database or look for the database if it already exists
    data_path_compressed = os.path.join(temp_repo_dir, mapping_info)  # Compressed .gz format
    data_path = os.path.splitext(data_path_compressed)[0]  # Not compressed
    if not os.access(data_path, os.R_OK) or not os.path.isfile(data_path):  # Check if it already exist
        subprocess.run((f"cd {temp_repo_dir}; curl -L -R -O ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/"
                        f"release_0.92/{mapping_info}"), shell=True)  # Download otherwise
        subprocess.run(f"cd {temp_repo_dir}; gzip -d {data_path_compressed}", shell=True)  # Uncompress the file
    data = dict()
    with open(data_path, "r") as input_handle:  # Open the file
        # Read the file, strip white spaces at the ends, split lines into cells
        mane_raw = [i.strip().split('\t') for i in input_handle.readlines() if not i.startswith('#')]
        for entry in mane_raw:  # For every entry in the data file
            if entry[9].startswith("MANE Select"):  # See if it is 'MANE Select' entry, ignore otherwise
                clean_entry = [i.split('.')[0] for i in [entry[1], entry[7]]]  # Remove version of gene and transcript
                assert clean_entry[0] not in data  # Make sure every gene is found once in the database
                data[clean_entry[0]] = clean_entry[1]  # Add to the new dictionary
    return data


def mane_gene_maps_with_defined_protein(ensembl_relase_object, temp_repo_dir):
    """
    This is to make sure all MANE Select transcripts of a gene coding a protein with a Uniprot ID.
    :param ensembl_relase_object: Created by ensembl_relase_object_creator() function
    :param temp_repo_dir: String. Absolute path of directory to save necessary files
    :return: Dictionary of Gene → Transcript.
    """
    # Get the MANE initiative's mapping and their selection
    gene_map_transcript_mane = mane_gene_transcript_map(temp_repo_dir)
    # Get the gene to Uniprot map to see if a gene has a Uniprot ID
    gene_uniprot = db_to_uniprot_mapping(temp_repo_dir, mapping="Ensembl")
    # Get Uniprot to transcript map
    uniprot_transcript = uniprot_to_db_mapping_multiple(temp_repo_dir, "Ensembl_TRS")
    # Get transcript to Uniprot map
    transcript_uniprot = db_to_uniprot_mapping(temp_repo_dir, mapping="Ensembl_TRS")
    data = dict()
    for gene in gene_map_transcript_mane:  # For each gene in the MANEs gene list
        if gene_map_transcript_mane[gene] in transcript_uniprot:  # Check if it has a Uniprot ID
            data[gene] = gene_map_transcript_mane[gene]  # Use the transcript
        else:  # The transcript does not has a Uniprot ID
            if gene in gene_uniprot:  # If gene does not have a Uniprot ID, then do nothing.
                valid_transcripts = list()  # Create list to use in following lines
                # Follow Gene → Uniprot → Transcripts to find the transcripts coding a protein with Uniprot ID
                transcript_ids = list(
                    set([m for t in gene_uniprot[gene] for m in uniprot_transcript[t]["Ensembl_TRS"]]))
                for transcript_id in transcript_ids:
                    try:
                        # Get an transcript object to investigate the transcript
                        transcript_object = ensembl_relase_object.transcript_by_id(transcript_id)
                        assert transcript_object.biotype == "protein_coding"  # Make sure transcript is protein coding
                        # If TSL (Transcript support level) is not available. 100 is as arbitrarily large number.
                        tsl = transcript_object.support_level if transcript_object.support_level else 100
                        # Find the length of the transcript
                        transcript_length = range_sum(transcript_object.coding_sequence_position_ranges)
                        # Add the valid transcript to the valid_transcripts list
                        valid_transcripts.append([tsl, transcript_length, transcript_id])
                    except (ValueError, AssertionError):
                        # ValueError: If ensembl_relase_object cannot find the transcript id
                        # AssertionError: If transcript is not protein coding
                        pass
                try:
                    # Sort first by TSL then sort by length. Select the transcript for this gene
                    data[gene] = sorted(valid_transcripts, key=lambda x: (x[0], -x[1]))[0][2]
                except IndexError:
                    # If valid_transcripts is empty.
                    pass
    return data


def download_gtf(temp_repo_dir):
    """
    Function is to download or find the file in temp_repo_dir.
    :param temp_repo_dir: String. Directory to download or find the file
    :return: String. Path of gtf or gtf.gz
    """
    gtf_file = "Homo_sapiens.GRCh38.100.gtf.gz"  # File name in the FTP server
    # Paths to download the database or look for the database if it already exists
    gtf_gz_path = os.path.join(temp_repo_dir, gtf_file)  # Compressed .gz format
    gtf_path = os.path.splitext(gtf_gz_path)[0]  # Not compressed
    if os.access(gtf_path, os.R_OK) or os.path.isfile(gtf_path):  # Check if the file exist
        return gtf_path  # Return the path if it exist
    elif os.access(gtf_gz_path, os.R_OK) or os.path.isfile(gtf_gz_path):
        return gtf_gz_path  # Return the compressed file if it exist
    else:  # Download otherwise
        subprocess.run((f"cd {temp_repo_dir}; "
                        f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-100/"
                        f"gtf/homo_sapiens/{gtf_file}"), shell=True)
        return gtf_gz_path  # Return the compressed file


def gtf_parser(gtf_path, verbose=False):
    """
    This is to convert gtf file into list of dictionary elements.
    :param verbose: Boolean. Report or not report the progress.
    :param gtf_path: String. Path of the gtf or gtf.gz file
    :return: List of dictionaries for each entry
    """
    # Column names for GTF file. Source: https://www.ensembl.org/info/website/upload/gff.html
    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame"]
    extension = os.path.splitext(gtf_path)[1]  # Check the extension to see if it is compressed or not
    if extension == ".gtf":  # If uncompressed
        with open(gtf_path, "r") as gtf_handle:  # Open and read the file
            gtf = gtf_handle.readlines()
    elif extension == ".gz":  # If compressed
        with gzip.open(gtf_path, "rt") as gtf_handle:  # Open with gzip library and read the file
            gtf = gtf_handle.readlines()
    else:  # Report unexpected file extension
        raise Exception(f"{bcolors.FAIL}GTF file does not have a proper extension.{bcolors.ENDC}")
    # Remove the titles, strip the trailing white spaces, split the line into cells
    gtf = [i.strip().split('\t') for i in gtf if not i.startswith('#')]
    for i in range(len(gtf)):  # Parse the entries
        # Parse the attributes column, column 9, of each line
        # noinspection PyTypeChecker
        entry = dict(  # Get the result as dictionary
            [[k, l.replace("\"", "")] for k, l in  # 4) Remove double quote from value
             [j.strip().split(" ", 1) for j in  # 3) For each attribute split key from value
              gtf[i][8].split(';')  # 1) Split first with ';'
              if j]])  # 2) See if element contains anything (last character is ';', to remove the artifact of it)
        entry.update(dict(zip(columns, gtf[i][:8])))  # Append dictionary with remaining information (Info in columns)
        # noinspection PyTypeChecker
        gtf[i] = entry  # Replace the value in the list
        if verbose and (i % 100 == 0 or i == len(gtf) - 1):  # Show the progress if 'verbose' is True
            progressBarForTerminal(i, len(gtf) - 1)
    return gtf


def gene_exon_ranges(ensembl_release_object, gene_id):
    """
    For a given gene, create a list of ranges for the exon boundaries
    :param ensembl_release_object: Created by ensembl_release_object_creator() function
    :param gene_id: String ensembl id
    :return: Nested sorted list
    """
    exons = ensembl_release_object.exon_ids_of_gene_id(gene_id)  # Get the exon IDs of the gene.
    # Note that some exons have CDS but some others do not have
    ranges = list()  # Initialize a list to fill up
    for exon in exons:  # for each exon of the gene's exons
        exon_object = ensembl_release_object.exon_by_id(exon)  # Create a exon object of pyensembl library
        ranges.append((exon_object.start, exon_object.end))  # Populate the list by start and end positions of exons
    return sorted(reduce_range_list(ranges))  # Sort, reduce, and return


def gene_cds_ranges(ensembl_release_object, gene_id):
    """
    For a given gene, create a list of ranges for CDS boundaries
    :param ensembl_release_object: Created by ensembl_release_object_creator() function
    :param gene_id: String ensembl id
    :return: Nested sorted list
    """
    transcripts = ensembl_release_object.transcript_ids_of_gene_id(gene_id)  # Get the transcripts IDs of the gene.
    ranges = list()  # Initialize a list to fill up
    for i in transcripts:  # For each transcript of a gene
        try:
            # Get the transcript CDS ranges
            transcript_cds = ensembl_release_object.transcript_by_id(i).coding_sequence_position_ranges
            # Raises ValueError for some transcripts which do not have a CDS
        except ValueError:
            continue  # Ignore these transcripts
        for cds_range in transcript_cds:  # For CDS ranges of a trancript
            ranges.append(cds_range)  # Append the gene CDS ranges
    # reduce_range_list is to remove the overlapping CDS ranges, and redundant ranges
    return sorted(reduce_range_list(ranges))  # Sort, reduce, and return


def gene_longest_cds_with_defined_protein(ensembl_release_object, gene_id):
    """
    For a given gene, find the transcripts with longest CDS, return the genomic ranges of it.
    Makes sure CDS ranges originate from transcripts, which are protein coding.
    :param ensembl_release_object: Created by ensembl_release_object_creator() function
    :param gene_id: String ensembl id
    :return: Nested sorted list
    """
    transcript_ids = ensembl_release_object.transcript_ids_of_gene_id(gene_id)  # Get the transcripts IDs of the gene.
    lengths = list()  # Initialize a list to fill up
    for i in transcript_ids:  # For each transcript of a gene
        try:
            # Get the transcript CDS ranges
            transcript_object = ensembl_release_object.transcript_by_id(i)
            transcript_cds = transcript_object.coding_sequence_position_ranges
            # Raises ValueError for some transcripts which do not have a CDS
        except ValueError:
            continue  # Ignore these transcripts
        if transcript_object.is_protein_coding:  # If transcript is protein coding
            lengths.append((range_sum(transcript_cds), i))  # Calculate CDS length, add (length, id) tuples to list
    # Sort the transcripts based on their lengths, chose the longest one, report the CDS ranges of it.
    return ensembl_release_object.transcript_by_id(sorted(lengths)[-1][1]).coding_sequence_position_ranges


def reduce_range_list(ranges):
    """
    Best way to return minimum number of ranges from a collection of ranges
    :param ranges: List of ranges (nested list), or list of tuples
    :return: List of lists
    """
    # Sorting the list based on the lower limit, and then based on size as tie breaker
    # Note: We need the bigger ranges to come before the smaller ranges
    ranges.sort(key=lambda pair: (pair[0], - pair[1]))
    reduced = []  # New set of ranges are stored here
    parent = ranges[0]  # Use a parent range to decide if a range is inside other ranges
    # This will for sure be part of the solution, because it is the largest, leftmost range
    reduced.append(ranges[0])  # Add it to the reduced list

    for x in ranges:  # for each entry in ranges
        if parent[0] <= x[0] and x[1] <= parent[1]:  # This range is completely within another range
            continue  # Ignore it
        elif x[0] <= parent[1]:  # This range is partially inside the parent range
            parent = [parent[0], x[1]]  # Set the parent to cover this two range
        else:  # This range is completely outside the parent range
            parent = x
        # If the range is completely or partially outside other ranges...
        # Place it here to avoid duplicate code
        reduced.append(x)

    return reduced  # Return the solution


def kai_mati_longest_transcripts_fetch(rscript, rdata_path, output_path):
    """
    This file is to fetch which transcripts which used in Mati-Kai's analysis.
    It runs a RScript in
    :param rscript:  In following directory: module_conservation/mati_kai_transcripts.R
    :param rdata_path:  It should be "All_sigmoidfits_HEK.RData".
    :param output_path:  String. Absolute path of resulting txt file
    :return:  List of transcript names
    """
    subprocess.run((
        f"{which('RScript')} {rscript} "  # Run RScript
        f"{rdata_path} "  # Using the RData object / file
        f"{output_path} "  # Save the file to this path
    ), shell=True)
    with open(output_path, "r") as input_handle:  # Open the file
        transcripts = [i.strip() for i in input_handle.readlines()]  # Read the content
    return list(set(transcripts))  # Get the unique transcripts and return them


def sync_dictionaries(dict1, dict2):
    """
    Syncronizes keys of two dictionaries. Keeps only keys which are found in both dictionaries
    :param dict1: Dictionary as dictionary 1
    :param dict2: Other dictionary as dictionary 2
    :return: Tuples of dictionaries
    """
    intersection = dict1.keys() & dict2.keys()  # Find the keys which are found in both dictionaries
    output_dict1 = {i: dict1[i] for i in intersection}  # Create a new dictionary from dict1 with only these keys
    output_dict2 = {i: dict2[i] for i in intersection}  # Create a new dictionary from dict2 with only these keys
    return output_dict1, output_dict2  # Return the filtered dictionaries


# End of the script
