"""

"""


import os
import re
import subprocess
import sys
import warnings
import psutil
import time
import itertools
# TODO: warnings.simplefilter('ignore', np.RankWarning)

import joblib
import h5py
from contextlib import closing
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing
from scipy.optimize import basinhopping
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from typing import Union
# from statsmodels.stats.proportion import proportion_confint  # todo
import pyensembl
import math
import pysam
from shutil import which, rmtree
from Bio.Seq import translate
from Bio import pairwise2
import matplotlib.pyplot as plt
import seaborn as sns
from functools import partial
from scipy.signal import find_peaks, peak_widths, peak_prominences


script_path_infrastructure = __file__


class Infrastructre:

    def __init__(self, temp_repo_dir, organism: str, exclude_genes: list = None, ensembl_release=102,
                 include_gene3d=False, verbose=True):

        self.temp_repo_dir = temp_repo_dir
        self.exclude_genes = [] if not exclude_genes else exclude_genes
        self.ensembl_release = ensembl_release
        self.organism = organism
        self.verbose = verbose

        self.script_path_infrastructure = script_path_infrastructure
        self.script_path_gene_info_db = os.path.abspath(os.path.join(os.path.dirname(self.script_path_infrastructure), "gene_info_database.R"))
        self.script_path_gene_info_names_db = os.path.abspath(os.path.join(os.path.dirname(self.script_path_infrastructure), "gene_info_names.R"))
        self.script_path_gene_info_uniprot_db = os.path.abspath(os.path.join(os.path.dirname(self.script_path_infrastructure), "gene_info_uniprot.R"))

        gene_info_database = biomart_mapping(self.temp_repo_dir, self.script_path_gene_info_db,
                                             self.ensembl_release, self.organism)
        gene_info_names = biomart_mapping(self.temp_repo_dir, self.script_path_gene_info_names_db,
                                          self.ensembl_release, self.organism)
        gene_info_uniprot = biomart_mapping(self.temp_repo_dir, self.script_path_gene_info_uniprot_db,
                                            self.ensembl_release, self.organism)

        gene_info_database = gene_info_database[gene_info_database["transcript_biotype"] == "protein_coding"]
        self.gene_list = sorted(np.unique(gene_info_database["ensembl_gene_id"].dropna()))
        self.gene_info = gene_class_dict_generate(self.temp_repo_dir, self.gene_list, gene_info_database,
                                                  gene_info_names, gene_info_uniprot, verbose=self.verbose)

        ero = ensembl_release_object_creator(self.temp_repo_dir, self.ensembl_release, self.organism)

        transcript_list = sorted(np.unique(gene_info_database["ensembl_transcript_id"].dropna()))
        self.protein_genome = ProteinGenome(self.temp_repo_dir, transcript_list, ero, verbose=self.verbose)

        if include_gene3d:
            self.script_path_gene3d_db = os.path.abspath(os.path.join(os.path.dirname(self.script_path_infrastructure), "gene3d.R"))
            self.gene3d_database = EnsemblDomain(self.temp_repo_dir, self.script_path_gene3d_db, self.protein_genome,
                                                 ero, ensembl_release=self.ensembl_release,
                                                 organism=self.organism, verbose=self.verbose)


    def create_gene_matrix(self, gene_id):
        pass
        # todo: integrates all information about the gene (n, l); n-feature, length
        # aa, nt,
        # traslatome_rpm, sixtymer_rpm..,
        # stalling sites, co-co sites, selective-ribosome,
        # uniprot, gene3d, conservation

    # annotation:
    # uniprot_annotations_ranges

    # conservation
    # conservation score -> ranges

    # get_riboseq, get_domain, get_conservation gibi function'lar yapılacakGE


def gene_class_dict_generate(temp_repo_dir, gene_list: list, gene_info_database: pd.DataFrame,
                             gene_info_names: pd.DataFrame, gene_info_uniprot: pd.DataFrame,
                             recalculate: bool = False, verbose: bool = True) -> dict:
    """
    This function creates a dictionary of Gene class elements.
    :param temp_repo_dir: Full path directory where temporary files will be stored.
    :param gene_list: List of gene IDs, for which to create the dictionary. Ideally, it should be sorted.
    :param gene_info_database: biomart_mapping() output for 'gene_info_database.R' script.
    :param gene_info_names: biomart_mapping() output for 'gene_info_names.R' script.
    :param gene_info_uniprot: biomart_mapping() output for 'gene_info_uniprot.R' script.
    :param recalculate: If True, it will calculate anyway.
    :param verbose: If True, it will print to stdout about the process computer currently calculates.
    :return: Dictionary of gene_id → Gene instance
    """

    # Get the absolute output file name.
    output_path = os.path.join(temp_repo_dir, "gene_info_database.joblib")
    # Check if there is already a calculated object saved before.
    if not os.access(output_path, os.R_OK) or not os.path.isfile(output_path) or recalculate:
        # If no file exists, or overwrite is 'True'.
        output = dict()  # Initiate an dictionary to fill up.
        print_verbose(verbose, f"Gene information dictionary are being created.")
        # For each gene in the gene list,
        for ind, gene_id in enumerate(gene_list):
            # Print out the current process to stdout as a progress bar.
            progress_bar(ind, len(gene_list) - 1, suffix=f"    {gene_id}", verbose=verbose)
            # Get the information from the pandas_data_frame objects for the gene of interest
            gi = gene_info_database[gene_info_database["ensembl_gene_id"] == gene_id]
            gi_names = gene_info_names[gene_info_names["ensembl_gene_id"] == gene_id]
            gi_uniprot = gene_info_uniprot[gene_info_uniprot["ensembl_gene_id"] == gene_id]
            # Use the information to create a Gene object using the Gene class.
            output[gene_id] = Gene(gi, gi_names, gi_uniprot)
        print_verbose(verbose, f"Results are being written to directory: {temp_repo_dir}")
        # Save the resulting filled dictionary into a joblib file to load without calculating again in next runs.
        joblib.dump(output, output_path)
        print_verbose(verbose, f"Done: {output_path}")
        return output  # Return the resulting dictionary.
    else:  # If the joblib file is found in the temp directory, or overwrite is 'False'.
        print_verbose(verbose, f"Gene information dictionary is found in path: {output_path}")
        return joblib.load(output_path)  # Load and return the resulting dictionary.


class Gene:
    """
    Stores all the relevant data for a given gene, prioritizes the transcripts.
    """

    def __init__(self, one_gene_df: pd.DataFrame, one_gene_names: pd.DataFrame, one_gene_uniprot: pd.DataFrame):
        """
        Init function. 'A=k[k["ensembl_gene_id"] == gene_id]' where k is biomart_mapping() output.
        See gene_class_dict_generate() for its implementation.
        :param one_gene_df: 'A', where biomart_mapping was run with gene_info_database.R
        :param one_gene_names: 'A', where biomart_mapping was run with gene_info_names.R
        :param one_gene_uniprot: 'A', where biomart_mapping was run with gene_info_uniprot.R
        """
        # Below the class variables are first set to a 'temp_x' because in the following assert statement they will
        # be tested whether there is only one unique value for this for a given gene
        temp_gene_id = one_gene_df["ensembl_gene_id"].unique()  # Get the gene ID
        self.gene_id = temp_gene_id[0]  # Assign to gene_id
        # Get and assign the gene names for the given gene_id. Get external_gene_name and external_synonym as one list.
        self.gene_names = np.unique(np.concatenate((one_gene_names["external_gene_name"].dropna(),
                                                    one_gene_names["external_synonym"].dropna())))
        # Get and assign the uniprot IDs for the given gene_id. Get uniprotswissprot and uniprotsptrembl as one list.
        self.uniprot_ids = np.unique(np.concatenate((one_gene_uniprot["uniprotswissprot"].dropna(),
                                                     one_gene_uniprot["uniprotsptrembl"].dropna())))
        temp_chromosome = np.unique(np.array(one_gene_df["chromosome_name"], dtype=str))  # Get the chromosome name
        self.chromosome = str(temp_chromosome[0])  # Assign to chromosome
        temp_strand = one_gene_df["strand"].unique()  # Get the strand
        self.strand = "-" if temp_strand[0] == -1 else "+"  # Assign to strand
        temp_start = one_gene_df["start_position"].unique()  # Get start position
        self.start = int(temp_start[0])  # Assign to end position
        temp_end = one_gene_df["end_position"].unique()  # Get start position
        self.end = int(temp_end[0])  # Assign to end position
        # As mentioned above, these elements should be unique for a given gene.
        assert all([len(i) == 1 for i in [temp_chromosome, temp_strand, temp_start, temp_end]])
        # Run prioritize transcripts and assign it to 'transcripts' namespace.
        self.transcripts = self.prioritize_transcripts(one_gene_df)

    @ staticmethod
    def prioritize_transcripts(gene_df: pd.DataFrame) -> pd.DataFrame:
        """
        This method aims to prioritize transcripts based on their functional relevance and amount of evidence
        supporting its presence.
        :param gene_df: In the same format as one_gene_df, which is described in 'init'.
        :return: Sorted data frame of transcripts with their scoring by different methods.
        """
        # Keep only the relevant information for transcripts.
        transcripts = gene_df.drop(["ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"],
                                   axis=1).drop_duplicates()
        # Note these columns comes from gene_info_database.R, so if you chance this 'R' script, change here as well.
        # As pd.Dataframe.sort_values is used, the sorting is based on alphabetical order. 'By chance', for all others
        # like mane, gencode, scores can be correctly sorted by just sorting alphabetically. However, for appris
        # 'P' scores are actually better than 'alt' scores so I just renamed.
        transcripts["transcript_appris"] = transcripts["transcript_appris"].replace(
            ["alternative1", "alternative2"], ["renamed_alternative1", "renamed_alternative2"])
        # Check there is only one MANE select for a gene; because if not, the sorting will fail.
        if "transcript_mane_select" in transcripts.columns:  # only human contains MANE annotation
            assert len(np.unique(transcripts["transcript_mane_select"].dropna())) <= 1
        # Based on my research, sorting the transcripts by the following order will always give the best possible
        # transcript in terms of functional relevance and amount of evidence supporting its presence.
        if "transcript_mane_select" in transcripts.columns:  # only human contains MANE annotation
            transcripts.sort_values(by=["transcript_mane_select", "transcript_appris", "transcript_gencode_basic",
                                        "transcript_tsl", "external_transcript_name"], inplace=True, ignore_index=True)
        else:
            transcripts.sort_values(by=["transcript_appris", "transcript_gencode_basic",
                                        "transcript_tsl", "external_transcript_name"], inplace=True, ignore_index=True)
        # See https://m.ensembl.org/info/genome/genebuild/transcript_quality_tags.html for more information
        return transcripts  # Return the sorted and narrowed dataframe


class ProteinGenome:
    """
    Store all transcripts sequence as well as corresponding proteins' sequence. This class is created to solve out some
    problems in the databases that cause my calculations to suffer. For example, for some transcripts that has
    undefined 3' or 5', the length of protein sequence does not match the length of transcript sequence. The methods in
    this class solves the problem by adding extra Ns to the transcript sequence or trims the nucleotide sequence; the
    reference point is always protein sequence. In addition to this, the class convert genomic coordinates to protein
    or transcript ranges and vice versa.
    """

    def __init__(self, temp_repo_dir: str, transcript_list: list, ensembl_release_object: pyensembl.Genome,
                 verbose: bool = True, recalculate: bool = False):
        """
        Init function.
        :param temp_repo_dir: Full path directory where temporary files will be stored.
        :param transcript_list: List of gene IDs, for which to create the dictionary. Ideally, it should be sorted.
        :param ensembl_release_object: Output of ensembl_release_object_creator()
        :param verbose: If True, it will print to stdout about the process computer currently calculates.
        :param recalculate: If True, it will calculate anyway.
        """
        self.transcript_list = sorted(transcript_list)  # Sort the transcript list in case it is not sorted and assign
        self.ensembl_release = ensembl_release_object.annotation_version  # Save the 'annotation_version'
        self.temp_repo_dir = temp_repo_dir  # Save the 'temp_repo_dir'
        self.verbose = verbose  # Save the whether 'verbose' activated or not
        self.recalculate = recalculate  # Save the whether 'recalculate' activated or not.
        # Get the absolute output file name.
        self.output_file_name = os.path.join(self.temp_repo_dir, "protein_genome_instance.joblib")

        try:
            # Check if there is already a calculated object saved before.
            assert os.access(self.output_file_name, os.R_OK) and os.path.isfile(self.output_file_name)
            if self.recalculate:  # If 'recalculate' is True,
                print_verbose(verbose, f"Saved file is found at the path but 'recalculate' is activated: "
                                       f"{self.output_file_name}.")
                raise AssertionError  # Raise the error to go to the except statement.
            # Load the saved content from the directory.
            loaded_content = load_joblib("ProteinGenome", self.output_file_name, self.verbose)
            # Check the current run is consistent with the saved file
            consistent_with = all([  # Make sure all below is True
                self.transcript_list == loaded_content.transcript_list,  # All transcripts are the same, and sorted.
                self.ensembl_release == loaded_content.ensembl_release,  # The same ensembl release is used.
            ])
            if not consistent_with:  # If there is a problem, the stored does not match with current run.
                print_verbose(verbose, f"There is at least one inconsistency between input parameters and "
                                       f"loaded content. Recalculating.")
                raise AssertionError  # Raise the error to go to the except statement.
            # Otherwise, just accept the database saved in previous run.
            self.db = loaded_content.db
        except (AssertionError, AttributeError, FileNotFoundError):  # If an error is raised.
            print_verbose(verbose, f"Protein genome mapping are being calculated.")
            self.db = self.calculate_transcript_mapping(ensembl_release_object)  # Calculate to get the mappings
            # Save the resulting filled dictionary into a joblib file to load without calculating again in next runs.
            save_joblib(self, self.output_file_name, self.verbose)

    @staticmethod
    def consistent_coding_ranges(transcript_object: pyensembl.Transcript) -> tuple:
        """
        Main function of the class purpose. Explanation is given as a part of class description above. The function
        makes the transcript and associated protein sequences consistent, makes number of nucleotides three times
        amino acid number. As a byproduct, this function also removes stop codon. Function outputs sequences as well as
        genomic position ranges in the same order (considering the strand) which facilitates the conversions between
        genome, protein and transcript positions.
        :param transcript_object: Created by ensembl_release_object_creator() function
        :return: Tuple of all relevant 'corrected' information for given transcript. Namely, (1) correctly ordered
        CDS genomic ranges without stop codon, (2) protein sequence, (3) CDS without stop codon, (4) chromosome name,
        (5) strand.
        """

        def protein_coding_ranges(tx_object: pyensembl.Transcript) -> list:
            """
            Some transcript has no defined start or stop codon. This causes the pyensembl module to fail in to find
            coding sequence with coding_sequence method. This function is to replace this method. It outputs a series
            of ranges without introns, and the positions are starting from 0. This positional ranges will be useful
            to get CDS sequence from spliced cDNA sequence of transcript.  Please note that
            below conditions are already tested. For a given transcript, all CDSs are at the same strand. For a given
            transcript, coding_sequence_position_ranges method gives CDSs in order. For transcripts which does not
            raise error with coding_sequence, this function gives identical results.
            :param tx_object: Created by ensembl_release_object_creator() function
            :return: Positions of coding sequence as series of ranges.
            """
            # Get the exons defined to the transcript as pyensembl.Exon objects
            transcript_object_exons = tx_object.exons
            # Make a database query to fetch exons which contains also CDS
            exon_numbers = tx_object.db.query(select_column_names=["exon_number"], feature="CDS",
                                              filter_column="transcript_id", filter_value=tx_object.id)
            first_coding_exon = int(exon_numbers[0][0]) - 1  # Which one is the first coding exon
            # Assign first coding exon to a variable
            origin_exon = transcript_object_exons[first_coding_exon]
            # Get the length of the transcript's exons until first coding exon.
            offset = ensembl_range_sum(tx_object.exon_intervals[:first_coding_exon])
            cdr = tx_object.coding_sequence_position_ranges  # Get CDS ranges for the transcript
            # Below if-else statement is to make CDS ranges starting from 0, so it subtracts the first CDS's genomic
            # starting position from all CDS ranges. However, introns are still there.
            if tx_object.strand == "+":
                cdr_relative_temp = [[s - origin_exon.start, e - origin_exon.start] for s, e in cdr]
            else:  # Notice it also switches start and end positions of the ranges
                cdr_relative_temp = [[origin_exon.end - e, origin_exon.end - s] for s, e in cdr]
            # Below while statement is to remove introns as well.
            cdr_relative = [cdr_relative_temp.pop(0)]  # Get the first item directly.
            while len(cdr_relative_temp) > 0:  # Keep until the last CDS ranges
                intron = cdr_relative_temp[0][0] - cdr_relative[-1][1] - 1  # Find intron length.
                # Subtract the intron length from all remaining CDS
                cdr_relative_temp = [[s - intron, e - intron] for s, e in cdr_relative_temp]
                # As the first item is corrected, get it.
                cdr_relative.append(cdr_relative_temp.pop(0))
            # Bring back skipped exons
            cdr_relative = [[offset + s, offset + e] for s, e in cdr_relative]
            # Coding sequence: "".join([transcript_object.sequence[s: e + 1] for s, e in cdr_relative])
            return cdr_relative

        def transcript_frame(transcript_as_query: str, protein_as_target: str, table: int) -> int:
            """
            Calculates the transcript frame by converting transcript (CDS of the transcript) into protein sequence
            then aligns to the protein. The frame with best alignment score will be returned at the end.
            :param transcript_as_query: Transcript (CDS of the transcript) nucleotide sequence
            :param protein_as_target: Protein amino acid sequence
            :param table: The conversion table for nucleotide to amino acid transformation. Integer number corresponds
            to the number at the following link: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi For human,
            standard table is '1', mitochondrial table is '2'.
            :return: Frame, either 1, 2 or 3.
            """
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')  # Translation and alignments throw warnings when success is so low.
                results = list()  # Initiate a list to keep alignment results.
                for frame_to_check in range(3):  # Test the alignment for frame 0, 1 and 2.
                    # Start transcript with appointed frame and then translate to protein sequence.
                    frame = translate(transcript_as_query[frame_to_check:], table=table)
                    # Try to align to the target protein sequence.
                    score = pairwise2.align.localxx(frame, protein_as_target, score_only=True)
                    if score:  # If a bit of alignment
                        results.append([score, frame_to_check])  # Append to the list.
                assert len(results) > 0, "No alignment!"  # Raise an error if there is no alignment at all.
                return sorted(results, reverse=True)[0][1]  # Sort based on alignment score and report the best frame.

        def find_overhang(query_seq: str, target_seq: str, reverse: bool, search_in: int = 10) -> int:
            """
            Compares query and target sequences' beginning (if reverse is False) or end (if reverse if True). It
            calculates how many nucleotides or amino acids more the target sequence is.
            :param query_seq: Nucleotides or amino acids sequence
            :param target_seq: Nucleotides or amino acids sequence
            :param reverse: To calculate 5' (False) or 3' end overhang
            :param search_in: Search space from beginning or end of sequences. The function is able to check the
            overhang with length half of search_in parameter. '10' is chosen arbitrarily after several rounds of
            trial and error, works well for Ensembl 102, Homo sapiens. Too big numbers causes some IndexError, too small
            numbers causes unsuccessful operation.
            :return: How much amino acids or nucleotides are overhung.
            """
            half_search = math.floor(search_in / 2)
            # Get first (or last if reverse is True) 'search_in' elements from both query and target.
            target_temp = target_seq[-search_in:][::-1] if reverse else target_seq[:search_in]
            query_temp = query_seq[-search_in:][::-1] if reverse else query_seq[:search_in]
            for query_offset in list(range(0, half_search)):
                # Get a slice from query with 'search_in / 2' length.
                query_last = query_temp[query_offset: query_offset + half_search]
                # Check if you can find the slice in the target (of first/last 'search_in' elements)
                relative_position = target_temp.find(query_last)
                if relative_position != -1:  # If the slice is found
                    return relative_position - query_offset  # Calculate and return number of the overhung elements.
                # Move to next iteration to get one element shifted slice.
            raise AssertionError("Error in find_overhang")  # If not overhung sequence is found, raise an error.

        try:  # For transcripts that are not coding, return relevant information here.
            # Get the CDS ranges of the transcript and convert it to a list.
            coding_position_ranges = [list(i) for i in transcript_object.coding_sequence_position_ranges]
        except ValueError:  # ValueError is raised when the transcript does not contain feature CDS.
            # Return relevant information, explained in the function's explanation above.
            # Exon ranges instead of CDS. No protein sequence. Whole exons sequence (whole cDNA) instead of CDS only.
            return ([[i, j] if transcript_object.strand == "+" else [j, i]
                     for i, j in transcript_object.exon_intervals],
                    None, transcript_object.sequence, transcript_object.contig, transcript_object.strand)

        if transcript_object.complete:  # If it has start and stop codons and a coding sequence divisible by 3.
            # Check below conditions anyway to make sure that the function returns correct result
            assert len(transcript_object.coding_sequence) == len(transcript_object.protein_sequence) * 3 + 3
            assert ensembl_range_sum(coding_position_ranges) == len(transcript_object.protein_sequence) * 3
            # Return relevant information, explained in the function's explanation above.
            return ([[i, j] if transcript_object.strand == "+" else [j, i] for i, j in coding_position_ranges],
                    transcript_object.protein_sequence, transcript_object.coding_sequence[:-3],
                    transcript_object.contig, transcript_object.strand)

        else:  # If the transcript is not complete, use above functions to correct the database information.
            target = transcript_object.protein_sequence  # Get the protein sequence
            query_positions = protein_coding_ranges(transcript_object)  # Get positions of CDS as list of ranges.
            # Get the coding sequence by using the above variable as explained before in protein_coding_ranges().
            query_transcript = "".join([transcript_object.sequence[s: e + 1] for s, e in query_positions])
            # Use a proper translation table according to the transcript's contig.
            translate_table = 2 if transcript_object.contig == "MT" else 1
            frame_start = transcript_frame(query_transcript, target, translate_table)  # Get the correct frame
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')  # Translation throws unnecessary warnings in some cases.
                # Get the protein sequence from transcripts CDS
                query = translate(query_transcript[frame_start:], table=translate_table)
            # Get rid of stop codon sign if only it is at the end.
            if query[-1] == "*":
                query = query[:-1]  # It is because amino acids like selenocysteine looks also '*'
            five_overhang = find_overhang(query, target, reverse=False)  # Calculate 5' overhang
            three_overhang = find_overhang(query, target, reverse=True)  # Calculate 3' overhang
            three_overhang_additional = ensembl_range_sum(coding_position_ranges) - frame_start - 3 * len(query)
            modify5 = -3 * five_overhang + frame_start  # Calculate how many nucleotides to add or remove
            modify3 = 3 * three_overhang - three_overhang_additional  # Calculate how many nucleotides to add or remove
            # Considering the transcript is at positive strand, change the first and end positions properly
            if transcript_object.strand == "+":
                coding_position_ranges[0][0] += modify5
                coding_position_ranges[-1][1] += modify3
            else:
                coding_position_ranges[0][1] -= modify5
                coding_position_ranges[-1][0] -= modify3
            # Add 'N' nucleotides to the overhung region, or trim out some nucleotides for both ends.
            query_transcript = query_transcript[-modify5:] if modify5 >= 0 else "N" * -modify5 + query_transcript
            query_transcript = query_transcript + "N" * modify3 if modify3 >= 0 else query_transcript[:modify3]
            # Check below conditions anyway to make sure that the function returns correct result
            assert ensembl_range_sum(coding_position_ranges) == len(transcript_object.protein_sequence) * 3
            assert len(query_transcript) == len(transcript_object.protein_sequence) * 3
            # Return relevant information, explained in the function's explanation above.
            return ([[i, j] if transcript_object.strand == "+" else [j, i] for i, j in coding_position_ranges],
                    transcript_object.protein_sequence, query_transcript,
                    transcript_object.contig, transcript_object.strand)

    def calculate_transcript_mapping(self, ensembl_release_object: pyensembl.Genome) -> dict:
        """
        Runs the consistent_coding_ranges() function for all transcripts in the transcript list of the object.
        :param ensembl_release_object: Output of ensembl_release_object_creator()
        :return: Transcript to corrected information dictionary.
        """
        # Make sure transcript list is unique
        assert len(np.unique(self.transcript_list)) == len(self.transcript_list)
        output = dict()  # Initiate a dictionary to fill up with transcript IDs as keys.
        for ind, transcript_id in enumerate(self.transcript_list):  # For each transcript in the list,
            # Print out the current process to stdout as a progress bar.
            progress_bar(ind, len(self.transcript_list) - 1, verbose=self.verbose)
            transcript_object = ensembl_release_object.transcript_by_id(transcript_id)  # Create a transcript object
            output[transcript_id] = self.consistent_coding_ranges(transcript_object)  # Call the function for the info
        return output

    def protein2genome(self, transcript_id: str, start: int, end: int) -> list:
        """

        :param transcript_id: Ensembl transcript ID without version.
        :param start: Start position of the protein range, query start position.
        :param end: End position of the protein range, query end position.
        :return: Nested list as genomic ranges
        """
        transcript = self.db[transcript_id]  # Get relevant transcript information from the database already prepared.
        # Check below condition to make sure the function will work successfully.
        assert 1 <= start <= end <= len(transcript[1]), f"Wrong range for transcript {transcript_id}: " \
                                                        f"min: 1, max: {len(transcript[1])}"
        # Convert amino acid positions to nucleotide positions. "-1"s are to make them python index.
        start_nt = start * 3 - 3
        end_nt = end * 3 - 1
        origin = transcript[0][0][0]  # The first nucleotide of the transcript.
        # Below if-else statement is to make CDS ranges starting from 0, so it subtracts the first CDS's genomic
        # starting position from all CDS ranges. Introns are also removed and saved.
        cdr_relative_temp = [[i - origin, j - origin] if transcript[4] == "+" else [origin - i, origin - j]
                             for i, j in transcript[0]]
        cdr_relative = [cdr_relative_temp.pop(0)]  # Get the first item directly.
        introns = [0]  # Initiate the list for intron lengths, first intron is assumed to be 0 for convenience.
        while len(cdr_relative_temp) > 0:  # Keep until the last CDS ranges
            intron = cdr_relative_temp[0][0] - cdr_relative[-1][1] - 1  # Find intron length.
            # Subtract the intron length from all remaining CDS
            cdr_relative_temp = [[s - intron, e - intron] for s, e in cdr_relative_temp]
            # As the first item is corrected, get it.
            cdr_relative.append(cdr_relative_temp.pop(0))
            introns.append(intron)  # Also save the intron length.
        # The logic is similar to what is done in protein_coding_ranges()
        introns = list(itertools.accumulate(introns))  # Calculate total introns until corresponding CDS.
        output_relative_introns = list()
        # Below is to calculate the genomic ranges, corresponding to the query protein range.
        for ind, (cdr_start, cdr_end) in enumerate(cdr_relative):
            if start_nt > cdr_end:
                continue  # Query range is not covered by this CDS range.
            elif cdr_start <= start_nt <= cdr_end and cdr_start <= end_nt <= cdr_end:
                output_relative_introns.append([start_nt + introns[ind], end_nt + introns[ind]])
                break
            elif cdr_start <= start_nt <= cdr_end < end_nt:
                output_relative_introns.append([start_nt + introns[ind], cdr_end + introns[ind]])
            elif start_nt <= cdr_start <= cdr_end < end_nt:
                output_relative_introns.append([cdr_start + introns[ind], cdr_end + introns[ind]])
            elif start_nt <= cdr_start <= end_nt <= cdr_end:
                output_relative_introns.append([cdr_start + introns[ind], end_nt + introns[ind]])
            else:
                break
        # Add the origin, the first nucleotide of the transcript, to the genomic ranges.
        # Output ranges are in the same order with the transcript.
        return [[i + origin, j + origin] if transcript[4] == "+" else [origin - i, origin - j]
                for i, j in output_relative_introns]


class EnsemblDomain:
    """
    This class has a data frame containing all protein domains of a given database. Each domain is appointed to
    a Ensembl protein ID. The class can return all domains of the database for a given transcript. The domain positions
    are kept as genomic ranges as well as protein positions. By means of ProteinGenome object, the conversion between
    protein, transcript and genome positions are quite easy.
    """

    # todo: The class can return all domains of a transcript

    def __init__(self, temp_repo_dir: str, rscript: str, protein_genome_instance: ProteinGenome,
                 ensembl_release_object: pyensembl.Genome, ensembl_release: int,
                 organism: str, verbose: bool = True, recalculate: bool = False):
        """
        Init function.
        :param temp_repo_dir: Full path directory where temporary files will be stored.
        :param rscript: Absolute path of R-Script, which fetches data from Ensembl. Check one of the script (e.g.
        gene3d.R) as a reference to properly use this class.
        :param protein_genome_instance: An instance of ProteinGenome class
        :param ensembl_release_object: Output of ensembl_release_object_creator()
        :param ensembl_release: Ensembl release version
        :param organism: Organism of interest to get the information from.
        :param verbose: If True, it will print to stdout about the process computer currently calculates.
        :param recalculate: If True, it will calculate anyway.
        """
        self.temp_repo_dir = temp_repo_dir
        self.rscript = rscript
        self.base_name = os.path.split(os.path.splitext(self.rscript)[0])[1]
        self.output_file_name = os.path.join(self.temp_repo_dir, f"{self.base_name}.joblib")
        self.ensembl_release = ensembl_release
        self.organism = organism
        self.verbose = verbose
        self.recalculate = recalculate

        try:
            # Check if there is already a calculated object saved before.
            assert os.access(self.output_file_name, os.R_OK) and os.path.isfile(self.output_file_name)
            if self.recalculate:  # If 'recalculate' is True,
                print_verbose(verbose, f"Saved file is found at the path but 'recalculate' is activated: "
                                       f"{self.output_file_name}.")
                raise AssertionError  # Raise the error to go to the except statement.
            # Load the saved content from the directory.
            loaded_content = load_joblib(f"EnsemblDomain for {self.base_name}", self.output_file_name, self.verbose)
            # Check the current run is consistent with the saved file
            consistent_with = all([  # Make sure all below is True
                self.base_name == loaded_content.base_name,  # The same domains are fetched.
                self.ensembl_release == loaded_content.ensembl_release  # The same ensembl release is used.
            ])
            if not consistent_with:  # If there is a problem, the stored does not match with current run.
                print_verbose(verbose, f"There is at least one inconsistency between input parameters and "
                                       f"loaded content. Recalculating.")
                raise AssertionError  # Raise the error to go to the except statement.
            # Otherwise, just accept the database saved in previous run.
            self.df = loaded_content.df
            self.columns = loaded_content.columns
        except (AssertionError, AttributeError, FileNotFoundError):  # If an error is raised.
            print_verbose(verbose, f"Ensembl domains are being calculated: {self.base_name}")
            self.df = biomart_mapping(self.temp_repo_dir, self.rscript,
                                      self.ensembl_release, self.organism)  # Get the annotations
            self.columns = self.df.columns
            self.convert2genome(protein_genome_instance, ensembl_release_object)  # Get genome coordinates as well.
            # Save the resulting data frame into a joblib file to load without calculating again in next runs.
            save_joblib(self, self.output_file_name, self.verbose)

    def convert2genome(self, protein_genome_instance: ProteinGenome, ensembl_release_object: pyensembl.Genome):
        """
        Calculates genomic coordinates of the domains.
        :param protein_genome_instance: An instance of ProteinGenome class
        :param ensembl_release_object: Output of ensembl_release_object_creator()
        :return: None. Function changes the object's main data frame.
        """
        # R-Script should be in such a format that satisfy below three lines.
        protein_ids = self.df[self.columns[0]]
        domain_starts = self.df[self.columns[1]]
        domain_ends = self.df[self.columns[2]]
        # Initiate lists to fill up with conversion calculations.
        coordinates_ranges, coordinates_contig = list(), list()
        # For each domain in the database
        for ind, (protein_id, domain_start, domain_end) in enumerate(zip(protein_ids, domain_starts, domain_ends)):
            # Print out the current process to stdout as a progress bar.
            progress_bar(ind, len(protein_ids) - 1, verbose=self.verbose)
            # Get the transcript ID corresponding to the protein ID
            transcript_id = ensembl_release_object.transcript_id_of_protein_id(protein_id)
            # Using the ProteinGenome object, calculate the genomic coordinates corresponding to the domains.
            coordinates_ranges.append(protein_genome_instance.protein2genome(transcript_id, domain_start, domain_end))
            coordinates_contig.append(protein_genome_instance.db[transcript_id][3])
        # Add these newly filled lists as a new columns to the object's main data frame.
        self.df["genome_chromosome"] = coordinates_contig
        self.df["genome_coordinate"] = coordinates_ranges


class RiboSeqAssignment:
    """
    This class calculates RiboSeq assignment for replicates from a list of SAM files. It also saves the assignment as
    joblib file to be used in later analysis. It can calculate RPKM and RPM for genes and gene positions.
    """

    def __init__(self, sam_paths: list, temp_repo_dir: str, riboseq_assign_at: Union[int, str], riboseq_assign_to: str,
                 riboseq_group: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict,
                 footprint_len: Union[int, list] = None, n_core: int = multiprocessing.cpu_count() - 2,
                 exclude_genes: list = None, verbose: bool = True, recalculate: bool = False):
        """
        Init function.
        :param sam_paths: List of absolute path of SAM files, which are replicates of each other.
        :param temp_repo_dir: Full path directory where temporary files will be stored.
        :param riboseq_assign_at: The position on the footprint where the footprint will be assigned at. '0' denotes to
        5' end 1st position and '1' denotes to 5' end 2nd position. '-1' corresponds to 3' 1st position, '-2'
        corresponds to 3' 2nd position. If it is "auto", then offsets will be calculated by the instance method
        called calculate_offset.
        :param riboseq_assign_to: There are 2 options for this, either "best_transcript" or "gene". One selects the
        best transcript, which is determined by Gene class for this gene. The latter selects all CDS of a given gene,
        sorts and merges them. Using "gene" is not recommended as other parts of the pipeline is prepared for
        "best_transcript" but not for "gene" for now.
        :param riboseq_group: Unique string to identify the RiboSeq assignment. e.g "cocoassembly_monosome"
        :param protein_genome_instance: An instance of ProteinGenome class
        :param n_core: Number of computer cores to be used.
        :param footprint_len: Keeps only the footprint with designated length. If None, keeps all footprint lengths.
        :param gene_info_dictionary: Output of gene_class_dict_generate() function
        :param verbose: If True, it will print to stdout about the process computer currently calculates.
        :param recalculate: If True, it will calculate anyway.
        """
        self.sam_paths = sorted(sam_paths)  # Sort the transcript list in case it is not sorted and assign
        self.temp_repo_dir = temp_repo_dir  # Save the 'temp_repo_dir'
        # Get the absolute output file name.
        self.output_file_name = os.path.join(temp_repo_dir, f"riboseq_{riboseq_group}.joblib")
        self.gene_list = sorted(gene_info_dictionary.keys())  # Sort the gene list in case it is not sorted and assign
        self.exclude_gene_list = list() if not exclude_genes else exclude_genes
        self.riboseq_assign_at = riboseq_assign_at
        self.riboseq_assign_to = riboseq_assign_to
        self.verbose = verbose  # Save the whether 'verbose' activated or not
        self.riboseq_group = riboseq_group  # Define riboseq_group variable for a function.
        self.recalculate = recalculate  # Save the whether 'recalculate' activated or not.
        self.footprint_len = footprint_len
        self.n_core = n_core
        self.f3, self.f5 = 33, 30  # Number of nucleotides in 3' and 5' ends of the genes

        try:
            # Check if there is already a calculated object saved before.
            assert os.access(self.output_file_name, os.R_OK) and os.path.isfile(self.output_file_name)
            if self.recalculate:  # If 'recalculate' is True,
                print_verbose(verbose, f"Saved file is found at the path but 'recalculate' is activated: "
                                       f"{self.output_file_name}.")
                raise AssertionError  # Raise the error to go to the except statement.
            # Load the saved content from the directory.
            loaded_content = load_joblib(f"RiboSeq assignment for {self.riboseq_group}",
                                         self.output_file_name, self.verbose)
            # Check the current run is consistent with the saved file
            consistent_with = all([  # Make sure all below is True
                all([os.path.basename(s) == os.path.basename(l)  # All SAM files have the same base name
                     for s, l in itertools.zip_longest(self.sam_paths, loaded_content.sam_paths)]),
                # Output file base name is the same. Also checks riboseq_group and riboseq_assign_to.
                os.path.basename(self.output_file_name) == os.path.basename(loaded_content.output_file_name),
                self.riboseq_assign_at == loaded_content.riboseq_assign_at,
                self.riboseq_assign_to == loaded_content.riboseq_assign_to,
                self.footprint_len == loaded_content.footprint_len,
                len(self.gene_list) == len(loaded_content.gene_list),  # There are the same amount of gene.
                all([g == l for g, l in zip(self.gene_list, loaded_content.gene_list)]),  # They are correctly sorted.
            ])
            if not consistent_with:  # If there is a problem, the stored does not match with current run.
                print_verbose(verbose, f"There is at least one inconsistency between input parameters and "
                                       f"loaded content. Recalculating.")
                raise AssertionError  # Raise the error to go to the except statement.
            # Otherwise, just accept the database saved in previous run.
            self.gene_assignments = loaded_content.gene_assignments
            self.flanking_assignments = loaded_content.flanking_assignments
            self.gene_lengths = loaded_content.gene_lengths
            self.total_assigned_gene = loaded_content.total_assigned_gene
            self.total_assigned = loaded_content.total_assigned
            self.length_distribution = loaded_content.length_distribution
            self.exclude_gene_list = loaded_content.exclude_gene_list
            self.footprint_len = loaded_content.footprint_len
            self.offset = loaded_content.offset

        except (AssertionError, AttributeError, FileNotFoundError):  # If an error is raised.

            # Get the mapping from chromosomes to list of genes that they carry.
            chromosome_gene = self.chromosomes_genes_matcher(gene_info_dictionary, self.gene_list)

            if self.riboseq_assign_at == "auto":  # If auto calculation of the offset is selected.
                self.gene_assignments, self.flanking_assignments, self.offset = self.calculate_offset(
                    self.n_core, chromosome_gene, protein_genome_instance, gene_info_dictionary,
                    footprint_len=self.footprint_len, abundance_threshold=0.01, verbose=self.verbose)
            elif isinstance(self.riboseq_assign_at, int):  # If offset is defined by the user
                self.gene_assignments, self.flanking_assignments, self.offset = self.fixed_offset(
                    chromosome_gene, protein_genome_instance, gene_info_dictionary,
                    footprint_len=self.footprint_len, verbose=self.verbose)
            else:  # Raise Error otherwise
                raise ValueError("Error in riboseq_assign_at")

            # Calculate gene lengths as well, which can be used in different applications.
            self.gene_lengths = {i: self.gene_assignments[i].shape[1] for i in self.gene_list}
            self.exclude_genes_calculate_stats(self.exclude_gene_list)  # Exclude genes, calculate RPM, RPKM etc.

            # Save the resulting object into a joblib file to load without calculating again in next runs.
            save_joblib(self, self.output_file_name, self.verbose)

    def exclude_genes_calculate_stats(self, exclude_gene_list: list):
        """
        Function removes all assignments and converts the values to NaN for genes in the gene list. Normally, when this
        method is first called, the exclude_gene_list is empty so the method just calculates the statistics. After
        saving the object, the function can be called by the user (or as a part of the script as seen in init method of
        RiboSeqExperiment class) to exclude the genes. Excluding some genes, due to a reason like over-expression etc,
        after the assignment is to make sure the object saved in the directory actually contains all assignments so
        that saved RiboSeqAssignment object can be used more flexibly.
        :param exclude_gene_list: List of Ensembl gene IDs.
        :return: Creates or updates total_assigned_gene, total_assigned and exclude_gene_list variables.
        """
        for gene_id in exclude_gene_list:  # For genes in the list
            empty_matrix = np.empty(self.gene_assignments[gene_id].shape)  # Create an empty array with the same shape
            empty_matrix.fill(np.nan)  # Fill the array with np.nan variables
            self.gene_assignments[gene_id] = empty_matrix  # Assign the empty array to the gene
            # Apply the same thing for self.gene_flanking
            empty_matrix = np.empty(self.flanking_assignments[gene_id].shape)
            empty_matrix.fill(np.nan)
            self.flanking_assignments[gene_id] = empty_matrix
            if gene_id not in self.exclude_gene_list:
                self.exclude_gene_list.append(exclude_gene_list)  # Save which genes were excluded in the class variable
        # Calculate the total number of assignments for each gene and for all. Mainly, to use in RPM, RPKM calculations
        self.total_assigned_gene = {i: np.nansum(self.gene_assignments[i], axis=1) for i in self.gene_list}
        self.total_assigned = np.nansum(np.array(list(self.total_assigned_gene.values())), axis=0)

    @staticmethod
    def calculate_length_distribution(footprint_asg_reps: list):
        """
        todo: here and the function
        :param footprint_asg_reps:
        :return:
        """
        result = list()
        for replicate in footprint_asg_reps:
            replicate_distribution = dict()
            for chromosome in replicate:
                lengths, counts = np.unique(replicate[chromosome][:, 1], return_counts=True)
                for ind, the_len in enumerate(lengths):
                    replicate_distribution[the_len] = replicate_distribution.get(the_len, 0) + counts[ind]
            result.append(replicate_distribution)
        return result

    @staticmethod
    def filter_by_footprint_abundance(length_distribution: list, abundance_threshold: float):
        """
        todo: here and the function
        :param length_distribution:
        :param abundance_threshold:
        :return:
        """
        total_assigned = [sum(i.values()) for i in length_distribution]  # Calculate total assignment for replicates
        # Calculate probability for each footprint length separately for each replicates, and filter based on that.
        pros = [{k: i[k] / total_assigned[ind] for k in i.keys()} for ind, i in enumerate(length_distribution)]
        pros_filtered = [{k: i[k] for k in i.keys() if i[k] > abundance_threshold} for ind, i in enumerate(pros)]
        return sorted({k for i in pros_filtered for k in i})

    def _calculate_offset_helper(self, footprint_asg_lst: list, chromosome_gene: dict,
                                 positions_gene: dict, positions_flanking: dict,
                                 protein_genome_instance: ProteinGenome, gene_info_dictionary: dict,
                                 total_assigned: int, f5_initial: int, f3_initial: int, figure_dir: str,
                                 footprint_len: int) -> tuple:
        """
        Used in by the auto-offset-calculating-method of this class so that the process can be completed in multiple
        cores at the same time for each footprint length separately.
        """
        # Initialize the resulting corrected dictionaries with zero values for all.
        result_k, result_l = dict(), dict()
        for gene_id in self.gene_list:
            result_k[gene_id] = np.zeros((len(footprint_asg_lst), len(positions_gene[gene_id])))
            result_l[gene_id] = np.zeros((len(footprint_asg_lst), self.f5 + self.f3))
            if gene_id in self.exclude_gene_list:
                # Create directly as a resulting object.
                result_k[gene_id][result_k[gene_id] == 0] = np.nan
                result_l[gene_id] = np.array([np.nan] * (self.f5 + self.f3))

        aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        fl = math.ceil(footprint_len / 3)  # Footprint length in amino acid
        min_rpkm = 1  # Ignore genes that has lower RPKM value of 1.
        palette = sns.color_palette("cubehelix", 3)
        plt.figure(figsize=(7, 3))
        offsets = dict()  # Initialize offset dictionary to fill up in the end

        # Below line is not tested. It is designed to work for sixtymers with kernel size of 3 and for monomers with
        # kernel size of 2. However, extensive testing is required.
        kernel = [1, 1, 1] if footprint_len > 41 else [1, 1]  # Warning: "Why 41" should be investigated further!
        first_position_ignore = max(-7, math.floor(fl / 3) + 1 - fl)
        last_position_ignore = -2 if first_position_ignore < -5 else -1
        igr_ends = [first_position_ignore, last_position_ignore]  # Check only relevant positions in the metafootprint profile

        for frm in [0, 1, 2]:  # For each frames for a given footprint length

            # Assign the length-frame combinations to CDSs and flanking regions
            a_cds = self.footprint_counts_to_genes(footprint_asg_lst, chromosome_gene, positions_gene,
                                                   footprint_len=footprint_len, frame={frm}, verbose=False)
            a_fla = self.footprint_counts_to_genes(footprint_asg_lst, chromosome_gene, positions_flanking,
                                                   footprint_len=footprint_len, frame={frm}, verbose=False)

            prot_expected = ""  # Initialize string to calculate expected amino acid probability
            counter = 0  # Initialize to see how many footprint is used to estimate the A-site
            # Initialize a matrix to fill with normalized amino acid amounts for each positions for each amino acids
            fill_matrix = np.zeros((fl, len(aa)))

            for ind, gene_id in enumerate(self.gene_list):
                if gene_id in self.exclude_gene_list:  # Skip if the gene is among excluded genes.
                    continue
                # Get the best transcript which is used to calculate positions_gene variable, at least 1 exists.
                best_transcript = gene_info_dictionary[gene_id].transcripts.iloc[0][0]
                protein_sequence = protein_genome_instance.db[best_transcript][1]  # Get the protein sequence
                # Get the sum intensities for all positions of each replicates
                dt = np.sum(a_cds[self.gene_list[ind]], axis=0)
                # Calculate RPKM and make sure it is above min_rpkm defined above.
                if np.sum(dt) / np.sum(total_assigned) * 1e9 / len(dt) > min_rpkm:
                    dt = dt / np.mean(dt)  # Normalization for positions takes place here
                    for nt_position, norm_count in enumerate(dt):  # Get position, intensity
                        # Calculate the amino acid position of the given nucleotide position.
                        aa_position = math.floor(nt_position / 3)
                        if norm_count == 0 or aa_position - fl < 0:
                            # If the normalized intensity is 0 or aa position is smaller than max footprint length
                            continue
                        else:
                            fpca = protein_sequence[aa_position - fl: aa_position]  # Get FootPrint Covered Area
                            assert len(fpca) == fl  # Make sure it is the same as expected footprint length
                            prot_expected += protein_sequence  # Save the protein sequence to get estimated probs.
                            counter += 1  # Increment the counter as this position will be useful for the calculation.
                            for i_fl in range(fl):  # For each position in footprint (as amino acids)
                                try:
                                    # Increase the relevant position with normalized intensity
                                    fill_matrix[i_fl, aa.index(fpca[i_fl])] += norm_count
                                except ValueError:
                                    pass  # Raises ValueError for non-conventional amino acids
            if counter < 1000:
                offsets[frm] = np.nan
            else:
                # Calculate the observed amino acid probabilities for each position
                dfk_observed = pd.DataFrame(fill_matrix / np.mean(np.sum(fill_matrix, axis=1)), columns=aa)
                # Calculate the expected amino acid probabilities for each position
                dfk_expected = pd.DataFrame([ProteinAnalysis(prot_expected).get_amino_acids_percent()] * fl)
                # Calculate KL divergence (entropy) for each position
                sonc = [stats.entropy(dfk_observed.iloc[i], dfk_expected.iloc[i]) for i in range(fl)]
                # Get the sum of consecutive positions (E-P-A positions) and find the maximum
                max_kl_couple = np.argmax(np.convolve(sonc, kernel, mode="valid")[igr_ends[0]: igr_ends[1]])
                a_site_pos = (igr_ends[0] + max_kl_couple) * 3  # Get the a-site position
                offsets[frm] = a_site_pos  # Add to the dictionary to report

                plt.plot(sonc, color=palette[frm], marker=".", label=frm)

                for gene_id in self.gene_list:  # For all gene in the gene list
                    if gene_id not in self.exclude_gene_list:  # Skip if the gene is among excluded genes
                        # Combine flanking and CDS together into one array
                        merged = np.hstack([a_fla[gene_id][:, :f5_initial],
                                            a_cds[gene_id],
                                            a_fla[gene_id][:, f5_initial:]])
                        # Get a slice considering final f5, f3 values and the offset.
                        merged = merged[:, f5_initial - self.f5 - a_site_pos: self.f3 - f3_initial - a_site_pos]
                        result_k[gene_id] += merged[:, self.f5: -self.f3]  # Save the result
                        result_l[gene_id] += np.hstack([merged[:, :self.f5], merged[:, -self.f3:]])  # Save the result

        plt.title(f"Length {footprint_len}")
        plt.legend(title='Frames', bbox_to_anchor=(1, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(figure_dir, f"{footprint_len}.pdf"))

        max_possible = 2 ** 16 - 1  # Maximum assignment to a single nucleotide position; as np.uint16 is used.
        for gene_id in self.gene_list:  # For all gene in the gene list
            if gene_id not in self.exclude_gene_list:  # Skip if the gene is among excluded genes
                assert np.max(result_k[gene_id]) <= max_possible
                assert np.max(result_l[gene_id]) <= max_possible
                result_k[gene_id] = result_k[gene_id].astype(np.uint16)
                result_l[gene_id] = result_l[gene_id].astype(np.uint16)

        return result_k, result_l, offsets

    def calculate_offset(self, n_core: int, chromosome_gene: dict, protein_genome_instance: ProteinGenome,
                         gene_info_dictionary: dict, footprint_len: Union[int, list],
                         abundance_threshold: float, verbose: bool) -> tuple:
        """
        # todo: comment here
        :param n_core:
        :param chromosome_gene:
        :param protein_genome_instance: An instance of ProteinGenome class
        :param gene_info_dictionary: Output of gene_class_dict_generate() function
        :param footprint_len:
        :param abundance_threshold:
        :param verbose:
        :return:
        """

        figure_dir = os.path.join(self.temp_repo_dir, f"riboseq_{self.riboseq_group}_figures_auto_offset")
        cleartarget(figure_dir)

        # Assuming each process takes 6 gb, calculate n_core to be used
        n_core = min(n_core, math.floor(psutil.virtual_memory().available / 6e+9))
        f5_initial, f3_initial = 120, 123  # Temporary values for f5 and f3 flanking positions.

        # Get array of genomic positions for all genes, create a dictionary out of it.
        positions_gene, positions_flanking = self.create_positions_genome(
            protein_genome_instance, gene_info_dictionary, f5_initial, f3_initial)

        # Assign footprints to genomic positions by calling footprint_assignment method for each SAM file.
        print_verbose(verbose, f"Footprints are being assigned to genomic coordinates.")
        footprint_asg_lst = [self.footprint_assignment(  # Assignment is at -1 position initially.
            sam_path, riboseq_assign_at=-1, split_by_length=True, calculate_distribution=False,
            footprint_len=footprint_len, verbose=verbose, only_dist=False)
            for sam_path in self.sam_paths]
        # Get the length distribution and save as the object variable
        self.length_distribution = self.calculate_length_distribution(footprint_asg_lst)
        # Calculate the abundant footprints to continue with
        selected_lengths = self.filter_by_footprint_abundance(self.length_distribution, abundance_threshold)
        selected_lengths = [i for i in selected_lengths if i >= 19]  # NOTE: Actually, depends on igr_ends
        total_assigned = [sum(list(i.values())) for i in self.length_distribution]  # To calculate RPKM in subprocess

        # Run individual gene assignments for these footprint lengths, assign at position -1,
        print_verbose(verbose, f"Gene assignments are being calculated for footprints "
                               f"with length {', '.join(map(str, selected_lengths))}.")
        ga_lol = partial(self._calculate_offset_helper, footprint_asg_lst, chromosome_gene,
                         positions_gene, positions_flanking, protein_genome_instance, gene_info_dictionary,
                         total_assigned, f5_initial, f3_initial, figure_dir)
        with closing(multiprocessing.Pool(n_core)) as executor:
            result_temp = executor.map(ga_lol, selected_lengths)

        # Save the offsets to report
        offsets = pd.DataFrame({i: r[2] for r, i in zip(result_temp, selected_lengths)}).T
        # Combine the resulting gene assignment results for flanking and CDS regions.
        gene_assignments, flanking_assignments = dict(), dict()
        max_possible = 2 ** 16 - 1  # Maximum assignment to a single nucleotide position; as np.uint16 is used.
        for gene_id in self.gene_list:  # For each gene in the gene list
            # Initiate the dictionary numpy array for flanking and CDS
            gene_assignments[gene_id] = np.zeros((len(footprint_asg_lst), len(positions_gene[gene_id])))
            flanking_assignments[gene_id] = np.zeros((len(footprint_asg_lst), self.f5 + self.f3))
            if gene_id in self.exclude_gene_list:  # If it is excluded gene, fill the np.nan values
                gene_assignments[gene_id][gene_assignments[gene_id] == 0] = np.nan
                flanking_assignments[gene_id][flanking_assignments[gene_id] == 0] = np.nan
            else:
                for r in result_temp:  # For each footprint length get the result
                    # Combine the results from all footprint lengths
                    gene_assignments[gene_id] += r[0][gene_id]
                    flanking_assignments[gene_id] += r[1][gene_id]
                assert np.max(gene_assignments[gene_id]) <= max_possible
                assert np.max(flanking_assignments[gene_id]) <= max_possible
                # Convert everything to np.uint16 so that it take up less space in the memory and the disk.
                gene_assignments[gene_id] = gene_assignments[gene_id].astype(np.uint16)
                flanking_assignments[gene_id] = flanking_assignments[gene_id].astype(np.uint16)

        print_verbose(verbose, f"The calculated offsets are as follows:\n{offsets}")
        return gene_assignments, flanking_assignments, offsets

    def fixed_offset(self, chromosome_gene: dict, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict,
                     footprint_len: Union[int, list], verbose: bool) -> tuple:
        """
        # todo: comment here
        :param chromosome_gene:
        :param protein_genome_instance: An instance of ProteinGenome class
        :param gene_info_dictionary: Output of gene_class_dict_generate() function
        :param footprint_len:
        :param verbose:
        :return:
        """

        # Get array of genomic positions for all genes, create a dictionary out of it.
        positions_gene, positions_flanking = self.create_positions_genome(
            protein_genome_instance, gene_info_dictionary, self.f5, self.f3)

        # Assign footprints to genomic positions by calling footprint_assignment method for each SAM file.
        print_verbose(verbose, f"Footprints are being assigned to genomic coordinates.")
        resulting_dicts = [self.footprint_assignment(
            sam_path, self.riboseq_assign_at, split_by_length=False, calculate_distribution=True,
            footprint_len=footprint_len, verbose=verbose, only_dist=False)
            for sam_path in self.sam_paths]
        footprint_genome_assignment_list, self.length_distribution = list(map(list, zip(*resulting_dicts)))

        # Assign footprints to gene positions by calling footprint_counts_to_genes method.
        print_verbose(verbose, f"Gene assignments are being calculated for footprints with all lengths.")
        gene_assignments = self.footprint_counts_to_genes(
            footprint_genome_assignment_list, chromosome_gene, positions_gene, verbose=False)
        flanking_assignments = self.footprint_counts_to_genes(
            footprint_genome_assignment_list, chromosome_gene, positions_flanking, verbose=False)
        offsets = None

        return gene_assignments, flanking_assignments, offsets

    def create_positions_genome(self, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict,
                                f5: int, f3: int) -> tuple:
        """
        Get array of genomic positions for each gene, create a dictionary out of it.
        :param protein_genome_instance: An instance of ProteinGenome class
        :param gene_info_dictionary: Output of gene_class_dict_generate() function
        :param f5: Number of nucleotides to take into account for 5' ends of the genes
        :param f3: Number of nucleotides to take into account for 3' ends of the genes
        :return: Tuple of dictionaries, which have gene_id as the key and genomic positions as the values. Note that
        the genomic positions already take into account the gene orientation. For example, if the gene in the reverse
        strand, the genomic positions will be decreasing numbers.
        """
        if self.riboseq_assign_to == "gene":
            # For all CDS of a given gene, sorts and merges them.
            # It is not recommended because almost all other methods, functions are written only for
            # the "best_transcript" for now.
            positions_gene = {gene_id: gene_entire_cds(protein_genome_instance, gene_info_dictionary, gene_id,
                                                       make_array=True) for gene_id in self.gene_list}
        elif self.riboseq_assign_to == "best_transcript":  # Only CDS of the best transcript, determined by Gene class
            positions_gene = {gene_id: best_transcript_cds(protein_genome_instance, gene_info_dictionary, gene_id,
                                                           make_array=True) for gene_id in self.gene_list}
        else:  # For now, there is only two methods for CDS selection.
            raise AssertionError("'riboseq_assign_to' variable must be either 'gene' or 'best_transcript'.")

        # Calculate positions_flanking to get the dictionary for flanking regions
        positions_flanking = dict()  # Initialize a dict
        for gene_id in positions_gene:  # Use the same genes which is used for gene assignment
            fp, lp = positions_gene[gene_id][0], positions_gene[gene_id][-1]  # Get first and last positions
            ff = 1 if fp < lp else -1  # Get the orientation of the genome positions
            near_start, near_end = np.arange(fp - f5 * ff, fp, ff), np.arange(lp + 1 * ff, lp + (f3 + 1) * ff, ff)
            positions_flanking[gene_id] = np.concatenate([near_start, near_end])
        # Resulting dictionary always follows 5'UTR, (no CDS), stop_codon,  3'UTR.

        return positions_gene, positions_flanking

    def metagene_at_ends(self, origin: str, window_length: int, min_rpkm: int,
                         gene_assignments: dict, flanking_assignments: dict, show_plot: bool) -> np.ndarray:
        """
        Calculates meta-gene profile for translation start or stop sites, which is particularly useful to see
        the characteristic periodicity of ribosome sequencing experiments.
        :param origin: Either 'start' or 'stop' to define the origin of the profile
        :param window_length: Indicates how much nucleotide is of interest
        :param min_rpkm: RPKM threshold to take genes into account for calculation of the profile.
        :param gene_assignments: Gene assignment dictionary, output of footprint_counts_to_genes
        :param flanking_assignments: Flanking assignment dictionary, output of footprint_counts_to_genes
        :param show_plot: Set True if the user wants to see the plot of the resuls as well
        :return: Average normalized read density for selected parameters.
        """
        dim = len(self.sam_paths)
        small_float = 1e-24

        # Calculate the total number of assignments for each gene and for all. To use in RPM, RPKM calculations
        total_assigned_gene = {i: np.nansum(gene_assignments[i], axis=1) for i in self.gene_list}
        total_assigned = np.nansum(np.array(list(total_assigned_gene.values())), axis=0).reshape(1, dim)

        # Create zero lists with the length of window_length
        the_sum, the_counts = np.zeros((dim, window_length)), np.zeros((dim, window_length))

        for gene_id in self.gene_list:  # For all genes in the gene_list
            if gene_id in self.exclude_gene_list:  # For excluded genes: np.nan values will distort the calculation
                continue  # Skip them
            rpm_pos_gene = (gene_assignments[gene_id] / total_assigned.T * 1e6)
            rpm_pos_flan = (flanking_assignments[gene_id] / total_assigned.T * 1e6)
            gene_rpkm = total_assigned_gene[gene_id] / total_assigned * 1e9 / rpm_pos_gene.shape[1]

            rpm_positions = np.hstack([rpm_pos_flan[:, :self.f5], rpm_pos_gene, rpm_pos_flan[:, -self.f3:]])  # Merge
            # If total RPM is zero, division for calculation of norm_rpm would be np.inf. As we already filter out
            # based on RPKM, RPM of the zero will be also filtered out. For here, zeros are replaced with a fairly small
            # number to suppress the warning message.
            total_rpm = np.nansum(rpm_positions, axis=1).reshape(dim, 1)
            total_rpm[total_rpm == 0] = small_float

            if origin == "start":
                norm_rpm = (rpm_positions / total_rpm)[:, :window_length]
                # For each replicate calculate if the replicate rpkm value for a gene is above threshold
                filter_rpkm = np.repeat((gene_rpkm > min_rpkm).reshape(dim, 1), norm_rpm.shape[1], axis=1).astype(int)
                # Calculates even though the gene length is less then window_length
                the_sum[:, :norm_rpm.shape[1]] += norm_rpm * filter_rpkm
                # Take into account the number of genes contributed to window positions. So, get the counts.
                the_counts[:, :norm_rpm.shape[1]] += np.ones((dim, norm_rpm.shape[1])) * filter_rpkm
            elif origin == "stop":  # The same procedure but for 'stop' origin
                norm_rpm = (rpm_positions / total_rpm)[:, -window_length:]
                filter_rpkm = np.repeat((gene_rpkm > min_rpkm).reshape(dim, 1), norm_rpm.shape[1], axis=1).astype(int)
                the_sum[:, -norm_rpm.shape[1]:] += norm_rpm * filter_rpkm
                the_counts[:, -norm_rpm.shape[1]:] += np.ones((dim, norm_rpm.shape[1])) * filter_rpkm
            else:
                raise ValueError
        metagene_profile = the_sum / the_counts

        if show_plot:  # If the user also want to see the profile visually.
            palette = sns.color_palette("cubehelix", len(self.sam_paths))
            plt.figure(figsize=(18, 4))
            mp_temp = metagene_profile * 1000
            for ind in range(len(self.sam_paths)):
                if origin == "start":
                    plt.plot(np.arange(-self.f5, len(mp_temp[ind]) - self.f5), mp_temp[ind],
                             color=palette[ind], marker=".", label=os.path.basename(self.sam_paths[ind]))
                else:  # origin == "stop":
                    plt.plot(np.arange(-len(mp_temp[ind]) + self.f3, self.f3), mp_temp[ind],
                             color=palette[ind], marker=".", label=os.path.basename(self.sam_paths[ind]))
            plt.title(f"Metagene profile around '{origin}' codon", fontweight="bold", y=1.05)
            plt.ylabel(r"Average normalized read density ($\times10^{-3}$)")
            plt.xlabel(f"Position from translation {origin} site")
            plt.legend(title='Replicates', bbox_to_anchor=(1, 1), loc='upper left')
            plt.tight_layout()
            plt.show()

        return metagene_profile  # Return the result as np.ndarray

    def plot_footprint_length_distribution(self, method="stack"):
        """
        Plots the footprint length distribution with following settings.
        :return: Always returns None
        """
        palette = sns.color_palette("cubehelix", len(self.sam_paths))
        the_lengths = sorted({lm for ld in self.length_distribution for lm in ld})
        if method == "stack":
            plt.figure(figsize=(12, 7.5))
            last_plotted = np.zeros(len(the_lengths))
            for ind, ld in enumerate(self.length_distribution):
                the_values = np.array([ld.get(fpl, 0) for fpl in the_lengths])
                plt.bar(the_lengths, the_values, bottom=last_plotted, color=palette[ind],
                        label=os.path.basename(self.sam_paths[ind]), linewidth=0)
                last_plotted += the_values
            plt.ylabel("Count")

        elif method == "percent":
            plt.figure(figsize=(12, 7.5))
            the_values = np.zeros((len(self.sam_paths), len(the_lengths)))
            for ind, ld in enumerate(self.length_distribution):
                the_values[ind] = np.array([ld.get(fpl, 0) for fpl in the_lengths])
            the_total = np.sum(the_values, axis=0)
            last_plotted = np.zeros(len(the_lengths))
            for ind, ld in enumerate(the_values):
                plt.bar(the_lengths, ld / the_total * 100, bottom=last_plotted, color=palette[ind],
                        label=os.path.basename(self.sam_paths[ind]), linewidth=0.5, edgecolor="white", width=1)
                last_plotted += ld / the_total * 100
            plt.ylabel("Percent")
        else:
            raise ValueError
        plt.xlabel("Footprint Length")
        plt.title("Genome-wide Footprint Length Distribution", fontweight="bold", y=1.05)
        plt.legend(title='Replicates', bbox_to_anchor=(1, 1), loc='upper left')
        plt.tight_layout()
        plt.gca().spines["top"].set_visible(False)
        plt.gca().spines["right"].set_visible(False)
        plt.show()

    @staticmethod
    def chromosomes_genes_matcher(gene_info_dictionary: dict, gene_list: list):
        """
        This script is to assign genes into chromosomes.
        :param gene_info_dictionary: Output of gene_class_dict_generate() function
        :param gene_list: List of genes to assign chromosomes
        :return: Dictionary of assigned genes, chromosome names as keys and the list of genes as values.
        """
        chromosomes = dict()  # Initialize a dictionary to fill up
        for gene_id in gene_list:  # For each gene in the gene_list
            contig = gene_info_dictionary[gene_id].chromosome  # Get the chromosome name
            # Find the chromosome name (contig), add the gene into the relevant list
            chromosomes[contig] = chromosomes.get(contig, []) + [gene_id]
        return chromosomes  # Return the resulting filled dictionary

    @staticmethod
    def footprint_assignment(sam_path: str, riboseq_assign_at: int, split_by_length: bool,
                             calculate_distribution: bool, only_dist: bool = False, verbose: bool = False,
                             footprint_len: Union[int, list] = None) -> Union[dict, tuple]:
        """
        Assigns the footprints into chromosomes and then into genomic positions. However, assigning to the genes are
        not a part of this method.
        :param verbose: If True, it will print to stdout about the process computer currently calculates.
        :param sam_path: UMI deduplicated final SAM file path to use for genome assignment.
        :param riboseq_assign_at: The position on the footprint where the footprint will be assigned at. '0' denotes to
        5' end 1st position and '1' denotes to 5' end 2nd position. '-1' corresponds to 3' 1st position, '-2'
        corresponds to 3' 2nd position. If it is a dict, then offsets could be different depending on footprint length.
        :param calculate_distribution: If True, calculated footprint length distributions will be calculated.
        :param split_by_length: If True, the assignment composes of two integers, one designating the position and the
        other shows the length of the footprint.
        :param only_dist: If True, only footprint length distribution will be calculated.
        :param footprint_len: Keeps only the footprint with designated length. If None, keeps all footprint lengths.
        :return: Dictionary of assigned positions as well as footprint length distribution
        """

        dist_single = dict()  # Initiate to calculate length distribution for given SAM
        assigned_positions = dict()  # Initialize dictionary with keys composed of chromosome names

        if verbose:  # Count the number of lines in SAM file to be able to print a progress bar.
            with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open SAM file with pysam library
                iteration = sam_handle.count()  # Call count method to find out total read count.

        counter = 0  # Initiate a counter to report the number of skipped footprints if exists.
        counter_footprint_length = 0  # To report the number of skipped footprints if footprint_len is not None.
        ache = "0123456789DIMNS"  # Allowed characters in cigar string.
        with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open SAM file with pysam library
            sam_iterator = sam_handle.fetch()  # Get the iterator to iterate through reads in for loop

            for ind, e in enumerate(sam_iterator):  # Iterate through the entries
                # Print out the current process to stdout as a progress bar.
                if verbose and (ind % 1000 == 0 or ind == iteration - 1):  # Print it in meaningful intervals
                    progress_bar(ind, iteration - 1, verbose=verbose)
                # This assertion is to make sure below if-else statement works perfectly well.
                assert all([i in ache for i in e.cigarstring]), "Cigar string must be composed of D, I, M, N, S only!"

                # Other than 'D', get_reference_positions method perfectly works. See the following link to understand
                # why we should keep 'D': https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/cigar.md
                if 'D' not in e.cigarstring:
                    # Gets the genomic coordinates of the aligned region
                    reference_positions = e.get_reference_positions()
                else:
                    # If there is a matching element (M) before to 'D', merge D to M.
                    for regex_match in re.finditer(r"([0-9]*)M[0-9]*D", e.cigarstring):
                        # regex_match.group() is like '48M2D'
                        m_len, d_len = re.findall(r"[0-9]+", regex_match.group())  # Get the numbers
                        m_len = int(m_len) + int(d_len)
                        e.cigarstring = e.cigarstring.replace(regex_match.group(), f"{m_len}M")
                    # Make sure that D's are merged to adjacent M's.
                    assert not re.search(r"([0-9]*D)", e.cigarstring),  "Make sure D's are merged to adjacent M"
                    reference_positions = e.get_reference_positions()  # Get genomic coordinates including 'D'

                # If only one footprint length is of interest.
                if footprint_len and isinstance(footprint_len, int) and len(reference_positions) != footprint_len:
                    counter_footprint_length += 1  # Count to report at the end.
                    continue  # Skip this footprint.
                # If a list of footprint length is of interest.
                if footprint_len and isinstance(footprint_len, list) and len(reference_positions) not in footprint_len:
                    counter_footprint_length += 1
                    continue

                # Fill up the dictionary for footprint length distribution
                if calculate_distribution:
                    dist_single[len(reference_positions)] = dist_single.get(len(reference_positions), 0) + 1
                if only_dist and calculate_distribution:  # Skip remaining if only distribution is of interest
                    continue

                # Calculate the assignment position for footprints in reverse or forward strand. Note that the
                # method normally accepts the assignment for forward strand.
                assignment_forward, assignment_reverse = riboseq_assign_at, -riboseq_assign_at - 1

                try:
                    if e.is_reverse:  # If reverse use the assigned position calculated in the beginning of the method.
                        assigned_nucleotide = reference_positions[assignment_reverse]
                    else:  # Please note that assignment_forward is actually equal to riboseq_assign_at variable.
                        assigned_nucleotide = reference_positions[assignment_forward]
                except IndexError:  # If list index out of range.
                    counter += 1  # This is the case if footprint is shorter than desired assigned position.
                    continue  # Skip this footprint.

                if split_by_length:
                    assigned_nucleotide = [assigned_nucleotide, len(reference_positions)]

                if e.reference_name not in assigned_positions:  # If the chromosome is not present in the dictionary
                    assigned_positions[e.reference_name] = [assigned_nucleotide]
                else:  # Add the assigned position of the footprint into the relevant list in the dictionary
                    assigned_positions[e.reference_name].append(assigned_nucleotide)

        if verbose and counter:  # Report the number of footprints skipped.
            print_verbose(verbose, f"{counter} footprints were skipped because assignment index was out of range.")

        if verbose and counter_footprint_length:  # Report the number of footprints skipped.
            print_verbose(verbose, f"{counter_footprint_length} footprints were skipped due to length constraint.")

        if only_dist and calculate_distribution:  # Return if only distribution is of interest
            return dist_single

        # After assignment ends for each entry in the chromosome dictionary
        for chromosome in assigned_positions:
            assigned_positions[chromosome] = np.array(assigned_positions[chromosome])  # Make np.array and sort

        if calculate_distribution:
            return assigned_positions, dist_single  # Return the filled dictionaries
        else:
            return assigned_positions

    def footprint_counts_to_genes(self, footprint_genome_assn_lst: list, chromosome_gene_map: dict,
                                  positions_map: dict, footprint_len: int = None, frame: set = None,
                                  verbose: bool = False) -> dict:
        """
        Assigns the footprints to genes. Output shows how many footprints are assigned at each position of a gene.
        :param footprint_genome_assn_lst: List of outputs of footprint_assignment function
        :param chromosome_gene_map: Output of chromosomes_genes_matcher function
        :param positions_map: Output of cds_ranges_to_array function
        :param verbose: If True, it will print to stdout about the process computer currently calculates.
        :param footprint_len: If not None, assign only the footprints with designated lengths
        :param frame: If not None, assign only the footprints with designated frame(s)
        :return: Dictionary mapping gene_id to numpy array (with the same length of gene), which compose of integers.
        """
        gene_footprint_assignment = dict()  # Initialize a dictionary to fill up

        for gene_id in self.gene_list:
            # Create an empty matrix to fill with raw counts, which originate from different replicates
            gene_footprint_assignment[gene_id] = np.zeros((len(footprint_genome_assn_lst), len(positions_map[gene_id])))

        # For each replicate do the same procedure
        for ind_assignment, footprint_genome_assignment in enumerate(footprint_genome_assn_lst):
            for ind, chr_name in enumerate(footprint_genome_assignment):  # For each chromosome
                # Print out the current process to stdout as a progress bar.
                progress_bar(ind, len(footprint_genome_assignment) - 1, verbose=verbose)
                # Get all genomic coordinates assigned to a footprint
                chr_pos_footprints = footprint_genome_assignment[chr_name]
                # Get the number of unique elements and the counts of these unique elements
                element_unique, element_count = np.unique(chr_pos_footprints, return_counts=True, axis=0)
                # Create a dictionary from this info, mapping genomic position with the number of footprints assigned.

                if not footprint_len:
                    # If no footprint_len is defined, chr_pos_footprints array must be in one dimensional, obtained by
                    # assignment method specifying split_by_length=False.
                    footprints_counts = dict(zip(element_unique, element_count))
                else:
                    # If footprint_len is defined, chr_pos_footprints array must be in two dimensional, obtained by
                    # assignment method specifying split_by_length=True.
                    temp_fc = np.hstack((element_unique, element_count.reshape(len(element_count), 1))).T
                    # Filter footprints by size, and discard the length row as they all the same
                    temp_fc = temp_fc[[0, 2], :][:, temp_fc[1] == footprint_len]
                    footprints_counts = dict(zip(temp_fc[0], temp_fc[1]))

                for gene_id in chromosome_gene_map.get(chr_name, []):  # For each gene in the same chromosome
                    pos_genes = positions_map[gene_id]  # Get the genomic positions of the gene
                    # For each position, look at the footprints_counts, get the number if found, 0 otherwise
                    temp_assignments = np.array([footprints_counts.get(i, 0) for i in pos_genes])
                    # Remove the count information from other frames, if frame is defined.
                    if frame is not None:
                        for frame_to_delete in list({0, 1, 2} - frame):
                            temp_assignments[frame_to_delete::3] = 0

                    # Make sure the line is not filled previously
                    assert np.max(gene_footprint_assignment[gene_id][ind_assignment]) == 0, "Error: Multiple genes"
                    gene_footprint_assignment[gene_id][ind_assignment] = temp_assignments  # Save the calculation

        max_possible = 2 ** 16 - 1  # Maximum assignment to a single nucleotide position; as np.uint16 is used.
        for gene_id in self.gene_list:
            if gene_id in self.exclude_gene_list:  # For genes in the list
                empty_matrix = np.empty(gene_footprint_assignment[gene_id].shape)  # An empty array with the same shape
                empty_matrix.fill(np.nan)  # Fill the array with np.nan variables
                gene_footprint_assignment[gene_id] = empty_matrix  # Assign the empty array to the gene
            else:
                assert np.max(gene_footprint_assignment[gene_id]) <= max_possible
                gene_footprint_assignment[gene_id] = gene_footprint_assignment[gene_id].astype(np.uint16)

        return gene_footprint_assignment  # Return filled dictionary

    def calculate_rpm_genes(self, gene_id: str, average: bool = True) -> Union[np.float64, np.ndarray]:
        """
        Calculates RPM for gene.
        :param gene_id: Ensembl gene ID, without version.
        :param average: To return the average of all replicates' calculation.
        :return: Calculated result in either np.float64 or np.ndarray, depending on 'average' parameter.
        """
        replicates = self.total_assigned_gene[gene_id] / self.total_assigned * 1e6
        return np.mean(replicates) if average else replicates

    def calculate_rpkm_genes(self, gene_id: str, average: bool = True) -> Union[np.float64, np.ndarray]:
        """
        Calculates RPKM for gene.
        :param gene_id: Ensembl gene ID, without version.
        :param average: To return the average of all replicates' calculation.
        :return: Calculated result in either np.float64 or np.ndarray, depending on 'average' parameter.
        """
        replicates = self.total_assigned_gene[gene_id] / self.total_assigned * 1e9 / self.gene_lengths[gene_id]
        return np.mean(replicates) if average else replicates

    def calculate_rpm_positions(self, gene_id: str, average: bool = True) -> np.ndarray:
        """
        Calculates RPM for each position in a gene.
        :param gene_id: Ensembl gene ID, without version.
        :param average: To return the average of all replicates' calculation.
        :return: RPMs of each positions. If average is True, then the array will have a shape of (gene_length, 1).
        Otherwise, the result has shape (gene_length, replicate_amount).
        """
        replicates = (self.gene_assignments[gene_id].T / self.total_assigned * 1e6).T
        return np.mean(replicates, axis=0) if average else replicates

    def export_h5(self, final_repo_dir: str, gene_info_dict: dict, verbose: bool = True):
        """
        Export the data as h5 file, which is required for some extarnal libraries; such as Ilia's R-Scripts
        :param final_repo_dir: The repository where the exported files should be saved.
        :param gene_info_dict: Output of gene_class_dict_generate() function.
        :param verbose: If True, it will print to stdout about the process computer currently calculates.
        """
        for sam_ind, sam_path in enumerate(self.sam_paths):  # Export each SAM file separately
            sam_name = os.path.basename(os.path.splitext(sam_path)[0])  # Get the base name for the next line
            output_file_name = os.path.join(final_repo_dir, f"riboseq_{self.riboseq_group}_{sam_name}_asite.h5")
            with h5py.File(output_file_name, "w") as data_file:  # Name output file properly, and create a h5 file
                for ind, gene_id in enumerate(self.gene_list):  # For each gene in the gene list
                    # Print out the current process to stdout as a progress bar.
                    prefix_bar = f"Exporting \'{sam_name}\' ({sam_ind+1}/{len(self.sam_paths)}): "
                    progress_bar(ind, len(self.gene_list) - 1, prefix=prefix_bar, verbose=verbose)
                    if gene_id in self.exclude_gene_list:
                        continue  # Skip if the gene is among excluded genes.
                    reduced = self.gene_assignments[gene_id][sam_ind]  # Get the assigned intensities for each position
                    positions = np.where(reduced != 0)  # Get the intensity positions where it is not 0
                    if len(positions[0]) == 0:
                        continue  # If there is not at least one assignment, skip the gene.
                    intensity = reduced[positions]  # Get relevant intensities, which corresponds to 'positions' array
                    arr = np.array([positions[0], intensity], dtype=np.uint32)  # Create an array out of those
                    if gene_info_dict[gene_id].strand == "-":
                        # In Ilia's method, the genes in the reverse strand was reversed. Just do the same.
                        arr = np.flip(arr, axis=1)
                    data_file[gene_id] = arr  # Write the result
                    # The metadata that are present also in Ilia's method
                    best_transcript_info = gene_info_dict[gene_id].transcripts.iloc[0].to_dict()
                    data_file[gene_id].attrs.update({
                        "cds_length": self.gene_lengths[gene_id],
                        "strand": gene_info_dict[gene_id].strand.encode('UTF-8'),
                        "chromosome": gene_info_dict[gene_id].chromosome.encode('UTF-8'),
                        "protein_id": best_transcript_info["ensembl_peptide_id"].encode('UTF-8'),
                        "transcript_id": best_transcript_info["ensembl_transcript_id"].encode('UTF-8'),
                        "transcript_name": best_transcript_info["external_transcript_name"].encode('UTF-8')})
                file_name = data_file.filename  # Record the name to report the output file path.
            print_verbose(verbose, f"Successfully exported: {file_name}")  # Report the output path


class RiboSeqExperiment:
    """
    Contains two RiboSeqAssignment, one is the experimental group and the other is for control (background)
    This is a parent class for many other experiment types such as CocoAssembly or Sixtymers, which contains their
    own specific methods.
    """

    def __init__(self, temp_repo_dir, sam_paths_background: list, sam_paths_experiment: list, name_experiment: str,
                 riboseq_assign_at: Union[int, str, list], riboseq_assign_to: str,
                 protein_genome_instance: ProteinGenome, gene_info_dictionary: dict,
                 footprint_len_experiment: Union[int, list] = None, footprint_len_background: Union[int, list] = None,
                 exclude_genes: list = None, n_core: int = multiprocessing.cpu_count() - 2,
                 verbose=True, recalculate=False):
        """
        Init function.
        :param temp_repo_dir: Full path directory where temporary files will be stored.
        :param sam_paths_background: Absolute paths of SAM files, which are replicates of each other, for background
        :param sam_paths_experiment: Absolute paths of SAM files, which are replicates of each other, for experiment
        :param name_experiment: Unique string to identify the RiboSeq experiment.  e.g "cocoassembly", "sixtymers"
        :param riboseq_assign_at: The position on the footprint where the footprint will be assigned at. '0' denotes to
        5' end 1st position and '1' denotes to 5' end 2nd position. '-1' corresponds to 3' 1st position, '-2'
        corresponds to 3' 2nd position. If it is "auto", then offsets will be calculated automatically. If a list of
        two elements are given, the first element in the list is for background and the second is for the experiment.
        :param riboseq_assign_to: There are 2 options for this, either "best_transcript" or "gene". One selects the
        best transcript, which is determined by Gene class for this gene. "gene" is not supported fully yet.
        :param protein_genome_instance: An instance of ProteinGenome class
        :param gene_info_dictionary: Output of gene_class_dict_generate() function
        :param footprint_len_experiment: Keeps only the footprint designated. If None, keeps all footprint lengths.
        :param footprint_len_background: Keeps only the footprint designated. If None, keeps all footprint lengths.
        :param exclude_genes: The Ensembl gene IDs that should be excluded from the calculations.
        :param n_core: Number of computer cores to be used.
        :param verbose: If True, it will print to stdout about the process computer currently calculates.
        :param recalculate: If True, it will calculate anyway.
        """
        self.temp_repo_dir = temp_repo_dir
        self.sam_paths_background = sam_paths_background
        self.sam_paths_experiment = sam_paths_experiment
        self.name_experiment = name_experiment
        self.riboseq_assign_to = riboseq_assign_to
        self.riboseq_assign_at = riboseq_assign_at if isinstance(riboseq_assign_at, list) \
            else [riboseq_assign_at, riboseq_assign_at]
        self.n_core = n_core
        self.exclude_genes = exclude_genes
        self.verbose = verbose
        self.recalculate = recalculate
        self.footprint_len_experiment = footprint_len_experiment
        self.footprint_len_background = footprint_len_background
        self.gene_list = sorted(gene_info_dictionary.keys())

        # Create the RiboSeqAssignment instances
        self.background = RiboSeqAssignment(self.sam_paths_background, self.temp_repo_dir, self.riboseq_assign_at[0],
                                            self.riboseq_assign_to, name_experiment + "_background",
                                            protein_genome_instance, gene_info_dictionary,
                                            footprint_len=self.footprint_len_background,
                                            exclude_genes=self.exclude_genes, n_core=self.n_core,
                                            recalculate=self.recalculate, verbose=self.verbose)
        self.experiment = RiboSeqAssignment(self.sam_paths_experiment, self.temp_repo_dir, self.riboseq_assign_at[1],
                                            self.riboseq_assign_to, name_experiment + "_experiment",
                                            protein_genome_instance, gene_info_dictionary,
                                            footprint_len=self.footprint_len_experiment,
                                            exclude_genes=self.exclude_genes, n_core=self.n_core,
                                            recalculate=self.recalculate, verbose=self.verbose)

        if self.exclude_genes:
            # After first instantiation of the class, the scripts saves the results as joblib files. Here, this chunk
            # of line was added because it allows users to exclude some more genes without recalculating (re-assigning)
            # from the beginning. Note that the saved RiboseqAssignment objects will not be changed after it is created.
            self.background.exclude_genes_calculate_stats(self.exclude_genes)
            self.experiment.exclude_genes_calculate_stats(self.exclude_genes)


class RiboSeqSixtymers(RiboSeqExperiment):
    """
    asd
    """

    def __init__(self, temp_repo_dir, sam_paths_background: list, sam_paths_experiment: list, name_experiment: str,
                 riboseq_assign_at: Union[int, str, list], riboseq_assign_to: str,
                 protein_genome_instance: ProteinGenome, gene_info_dictionary: dict,
                 footprint_len_experiment: Union[int, list] = None, footprint_len_background: Union[int, list] = None,
                 exclude_genes: list = None, n_core: int = multiprocessing.cpu_count() - 2,
                 verbose=True, recalculate=False):

        super().__init__(temp_repo_dir=temp_repo_dir,
                         sam_paths_background=sam_paths_background,
                         sam_paths_experiment=sam_paths_experiment,
                         name_experiment=name_experiment,
                         riboseq_assign_at=riboseq_assign_at,
                         riboseq_assign_to=riboseq_assign_to,
                         protein_genome_instance=protein_genome_instance,
                         gene_info_dictionary=gene_info_dictionary,
                         footprint_len_experiment=footprint_len_experiment,
                         footprint_len_background=footprint_len_background,
                         exclude_genes=exclude_genes,
                         n_core=n_core,
                         verbose=verbose,
                         recalculate=recalculate)

    def stalling_peaks_arpat(self, gene_id, mmc_threshold=1, normalized_peak_count_thr=5, get_top=5):
        # Arbitrarily 5 for all from the Arpat paper.
        land = self.experiment.calculate_rpm_positions(gene_id, average=True)  # library size normalized disome peaks

        normalized_peak_count = np.sum(land > 0)  # to test whether normalized peak count > 5
        mmc = self.background.calculate_rpm_genes(gene_id, average=True)  # mean monosome count
        if mmc < mmc_threshold or normalized_peak_count < normalized_peak_count_thr:
            return np.array([])
        else:
            normalized_rpm = land / mmc
            return normalized_rpm.argsort()[-get_top:][::-1]  # Arbitrarily 5

    def stalling_peaks_inecik_1(self, gene_id, percentile=90, window="hanning", window_len=23,
                                min_rpkm_sixtymers=-1, min_rpkm_background=1):
        try:
            assert type(self.aux_var_inecik_1)
            assert self.aux_var_inecik_1_params == (window, window_len, min_rpkm_sixtymers, min_rpkm_background)
        except (AttributeError, AssertionError):
            self.aux_var_inecik_1_params = (window, window_len, min_rpkm_sixtymers, min_rpkm_background)
            self.aux_var_inecik_1 = [exp_rpm_s[find_peaks(exp_rpm_s)[0]] for exp_rpm_s in
                                     self.inecik_gene_iterator(window, window_len, min_rpkm_sixtymers, min_rpkm_background)]
            self.aux_var_inecik_1 = np.log10(np.concatenate(self.aux_var_inecik_1, axis=0))
        threshold = np.percentile(self.aux_var_inecik_1, percentile)

        try:
            exp_rpm_s = self.inecik_get_gene(gene_id, window, window_len, min_rpkm_sixtymers, min_rpkm_background)
            peaks = find_peaks(exp_rpm_s)[0]
            peak_values = exp_rpm_s[peaks]
            return peaks[np.log10(peak_values) > threshold] - 1  # Added zero in the beginning
        except AssertionError:
            return np.array([])

    def stalling_peaks_inecik_2(self, gene_id, percentile=90, window="hanning", window_len=23, wlen=300,
                                min_rpkm_sixtymers=-1, min_rpkm_background=1):
        try:
            assert type(self.aux_var_inecik_2)
            assert self.aux_var_inecik_2_params == (window, window_len, wlen, min_rpkm_sixtymers, min_rpkm_background)
        except (AttributeError, AssertionError):
            self.aux_var_inecik_2_params = (window, window_len, wlen, min_rpkm_sixtymers, min_rpkm_background)
            self.aux_var_inecik_2 = list()
            for exp_rpm_s in self.inecik_gene_iterator(window, window_len, min_rpkm_sixtymers, min_rpkm_background):
                peaks, _ = find_peaks(exp_rpm_s)
                calc_prominences = peak_prominences(exp_rpm_s, peaks=peaks, wlen=wlen)
                calc_widths = peak_widths(exp_rpm_s, rel_height=1, peaks=peaks, prominence_data=calc_prominences, wlen=wlen)
                self.aux_var_inecik_2.append(calc_prominences[0] * calc_widths[0])
            self.aux_var_inecik_2 = np.log10(np.concatenate(self.aux_var_inecik_2, axis=0))
        threshold = np.percentile(self.aux_var_inecik_2, percentile)

        try:
            exp_rpm_s = self.inecik_get_gene(gene_id, window, window_len, min_rpkm_sixtymers, min_rpkm_background)
            peaks, _ = find_peaks(exp_rpm_s)
            calc_prominences = peak_prominences(exp_rpm_s, peaks=peaks, wlen=wlen)
            calc_widths = peak_widths(exp_rpm_s, rel_height=1, peaks=peaks, prominence_data=calc_prominences, wlen=wlen)
            calc = calc_prominences[0] * calc_widths[0]
            return peaks[np.log10(calc) > threshold] - 1  # Added zero in the beginning
        except AssertionError:
            return np.array([])

    def stalling_peaks_inecik_3(self, gene_id, probability=0.02, window="hanning", window_len=23, wlen=300,
                                min_rpkm_sixtymers=-1, min_rpkm_background=1):

        try:
            assert type(self.aux_var_inecik_3)
            assert self.aux_var_inecik_3_params == (window, window_len, wlen, min_rpkm_sixtymers, min_rpkm_background)
        except (AttributeError, AssertionError):
            self.aux_var_inecik_3_params = (window, window_len, wlen, min_rpkm_sixtymers, min_rpkm_background)
            calc_peak_width, calc_peak_prominence = list(), list()
            for exp_rpm_s in self.inecik_gene_iterator(window, window_len, min_rpkm_sixtymers, min_rpkm_background):
                peaks, _ = find_peaks(exp_rpm_s)
                calc_prominences = peak_prominences(exp_rpm_s, peaks=peaks, wlen=wlen)
                calc_widths = peak_widths(exp_rpm_s, rel_height=1, peaks=peaks, prominence_data=calc_prominences, wlen=wlen)
                calc_peak_prominence.extend(list(calc_prominences[0]))
                calc_peak_width.extend(list(calc_widths[0]))

            calc_peak_width, calc_peak_prominence = np.array(calc_peak_width), np.array(calc_peak_prominence)
            calc_filter_out = np.log(calc_peak_prominence) > -18
            calc_peak_width = calc_peak_width[calc_filter_out]
            calc_peak_prominence = calc_peak_prominence[calc_filter_out]

            params_h_additional = calc_peak_prominence.max()
            log_calc_peak_prominence = -np.log(calc_peak_prominence / calc_peak_prominence.max())
            params_h = stats.fisk.fit(log_calc_peak_prominence)

            filt_calc_peak_width = calc_peak_width[calc_peak_width != window_len - 1]
            params_w = stats.fisk.fit(filt_calc_peak_width)

            self.aux_var_inecik_3 = (params_w, params_h, params_h_additional)

        try:
            exp_rpm_s = self.inecik_get_gene(gene_id, window, window_len, min_rpkm_sixtymers, min_rpkm_background)
            peaks, _ = find_peaks(exp_rpm_s)
            calc_prominences = peak_prominences(exp_rpm_s, peaks=peaks, wlen=wlen)
            calc_widths = peak_widths(exp_rpm_s, rel_height=1, peaks=peaks, prominence_data=calc_prominences, wlen=wlen)
            # probability_prominence
            log_calc_peak_prominence = -np.log(calc_prominences[0] / self.aux_var_inecik_3[2])
            peak_prominence_probs = stats.fisk.cdf(log_calc_peak_prominence, *self.aux_var_inecik_3[1])  # because smaller values are higher peaks
            # probability_width
            peak_witdhs_probs = 1 - stats.fisk.cdf(calc_widths[0], *self.aux_var_inecik_3[0])
            bivariate_cumulative = peak_prominence_probs * peak_witdhs_probs
            return peaks[bivariate_cumulative < probability] - 1  # Added zero in the beginning
        except AssertionError:
            return np.array([])

    def stalling_peaks_inecik_4(self, gene_id, percentile=90, window="hanning", window_len=23, wlen=300,
                                min_rpkm_sixtymers=-1, min_rpkm_background=1):
        try:
            assert type(self.aux_var_inecik_4)
            assert self.aux_var_inecik_4_params == (window, window_len, wlen, min_rpkm_sixtymers, min_rpkm_background)
        except (AttributeError, AssertionError):
            self.aux_var_inecik_4_params = (window, window_len, wlen, min_rpkm_sixtymers, min_rpkm_background)
            self.aux_var_inecik_4 = list()
            for exp_rpm_s in self.inecik_gene_iterator(window, window_len, min_rpkm_sixtymers, min_rpkm_background):
                peaks, _ = find_peaks(exp_rpm_s)
                calc_prominences = peak_prominences(exp_rpm_s, peaks=peaks, wlen=wlen)
                self.aux_var_inecik_4.extend(list(calc_prominences[0]))
            self.aux_var_inecik_4 = np.log10(np.array(self.aux_var_inecik_4))
        threshold = np.percentile(self.aux_var_inecik_4, percentile)

        try:
            exp_rpm_s = self.inecik_get_gene(gene_id, window, window_len, min_rpkm_sixtymers, min_rpkm_background)
            peaks, _ = find_peaks(exp_rpm_s)
            calc_prominences = peak_prominences(exp_rpm_s, peaks=peaks, wlen=wlen)
            return peaks[np.log10(calc_prominences[0]) > threshold] - 1  # Added zero in the beginning
        except AssertionError:
            return np.array([])

    def stalling_peaks_inecik_5(self, gene_id, probability=0.02, min_rpkm_sixtymers=-1, min_rpkm_background=1):

        try:
            assert type(self.aux_var_inecik_5)
            assert self.aux_var_inecik_5_params == (min_rpkm_sixtymers, min_rpkm_background)
        except (AttributeError, AssertionError):
            self.aux_var_inecik_5_params = (min_rpkm_sixtymers, min_rpkm_background)
            aux_var_temp = list()
            for exp_rpm_s in self.inecik_gene_iterator(None, None, min_rpkm_sixtymers, min_rpkm_background,
                                                       smoothen=False):
                aux_var_temp.extend(list(exp_rpm_s))
            aux_var_temp = np.array(aux_var_temp)
            aux_var_temp = aux_var_temp[aux_var_temp != 0]

            params_additional = aux_var_temp.max()
            aux_var_temp = -np.log(aux_var_temp / params_additional)
            params = stats.fisk.fit(aux_var_temp)
            #x_data_mock = np.linspace(aux_var_temp.min(), aux_var_temp.max(), 512)
            #y_data_fit = stats.fisk.pdf(x_data_mock, *params)
            #plt.hist(aux_var_temp, bins=512, density=True, alpha=0.5)
            #plt.plot(x_data_mock, y_data_fit, color="black", linestyle='--')
            #plt.show()
            del aux_var_temp
            self.aux_var_inecik_5 = (params, params_additional)
        try:
            exp_rpm_s = self.inecik_get_gene(gene_id, None, None, min_rpkm_sixtymers, min_rpkm_background, smoothen=False)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                log_calc = -np.log(exp_rpm_s / self.aux_var_inecik_5[1])
                log_calc[log_calc == np.inf] = np.nan
            props = stats.fisk.cdf(log_calc, *self.aux_var_inecik_5[0])  # because smaller values are higher peaks
            return np.where(props < probability)[0]
        except AssertionError:
            return np.array([])

    #   todo: her bir replicate için ayrı hesaplamalı!
    #   todo: sonra logo bulma işine gir

    def inecik_get_gene(self, gene_id, window, window_len, min_rpkm_sixtymers, min_rpkm_background, smoothen=True):
        exp_rpkm = self.experiment.calculate_rpkm_genes(gene_id)
        bac_rpkm = self.background.calculate_rpkm_genes(gene_id)
        exp_rpm_bs = self.experiment.calculate_rpm_positions(gene_id)
        assert exp_rpkm > min_rpkm_sixtymers and bac_rpkm > min_rpkm_background
        rpm_gene = self.background.calculate_rpm_genes(gene_id)
        if smoothen:
            assert len(exp_rpm_bs) > window_len
            exp_rpm_s = [0] + list(smooth_array(exp_rpm_bs, window_len=window_len, window=window) / rpm_gene) + [0]
            return np.array(exp_rpm_s)
        else:
            exp_rpm_s = exp_rpm_bs / rpm_gene
            return exp_rpm_s

    def inecik_gene_iterator(self, window, window_len, min_rpkm_sixtymers, min_rpkm_background, smoothen=True):
        for gene_id in self.gene_list:
            try:
                yield self.inecik_get_gene(gene_id, window, window_len, min_rpkm_sixtymers,
                                           min_rpkm_background, smoothen=smoothen)
            except AssertionError:
                continue

    def plot_result(self, smoothen, gene_id, function, *args, **kwargs):
        total_exp = self.experiment.calculate_rpkm_genes(gene_id)
        total_tra = self.background.calculate_rpkm_genes(gene_id)
        rpm_tra = self.background.calculate_rpm_genes(gene_id)
        if smoothen:
            arr = smooth_array(self.experiment.calculate_rpm_positions(gene_id), window_len=23, window="hanning") / rpm_tra
        else:
            arr = self.experiment.calculate_rpm_positions(gene_id) / rpm_tra
        peaks = function(gene_id, *args, **kwargs)
        # Plot
        fig, ax = plt.subplots(1, 1, figsize=(14, 4))
        fig.suptitle(gene_id, y=1.1, fontweight="bold")
        ax.plot(arr, alpha=1, color="salmon")
        ax.scatter(peaks, [arr[p] for p in peaks], color="black", alpha=1, s=25)
        ax.set_ylim(0, arr.max() * 1.35)
        ax.text(0.01, 0.99, f"{gene_id} / Index {self.gene_list.index(gene_id)}\n"
                            f"RPKM Sixtymers: {round(total_exp, 2)}\n"
                            f"RPKM Translatome: {round(total_tra, 2)}",
                      fontsize=6, transform=ax.transAxes, verticalalignment='top', horizontalalignment="left")
        ax.text(0.99, 0.99, f"Max: {round(arr.max(), 5)}\nPeaks: {len(peaks[0])}", fontsize=6,
                      transform=ax.transAxes, verticalalignment='top', horizontalalignment="right")
        ax.axes.get_yaxis().set_visible(False)
        ax.tick_params(labelsize=6)
        return arr, peaks


class RiboSeqSelective(RiboSeqExperiment):

    def __init__(self, temp_repo_dir, sam_paths_background: list, sam_paths_experiment: list, name_experiment: str,
                 riboseq_assign_at: Union[int, str, list],
                 riboseq_assign_to: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict,
                 footprint_len_experiment: Union[int, list] = None, footprint_len_background: Union[int, list] = None,
                 exclude_genes: list = None, n_core: int = multiprocessing.cpu_count() - 2,
                 verbose=True, recalculate=False):
        super().__init__(temp_repo_dir=temp_repo_dir,
                         sam_paths_background=sam_paths_background,
                         sam_paths_experiment=sam_paths_experiment,
                         name_experiment=name_experiment,
                         riboseq_assign_at=riboseq_assign_at,
                         riboseq_assign_to=riboseq_assign_to,
                         protein_genome_instance=protein_genome_instance,
                         gene_info_dictionary=gene_info_dictionary,
                         footprint_len_experiment=footprint_len_experiment,
                         footprint_len_background=footprint_len_background,
                         exclude_genes=exclude_genes,
                         n_core=n_core,
                         verbose=verbose,
                         recalculate=recalculate)

    def calculate_binding_positions(self, normalized_rpm_threshold, min_gene_rpm_background, min_gene_rpkm_background):
        pass


class RiboSeqCoco(RiboSeqExperiment):

    def __init__(self, temp_repo_dir, sam_paths_monosome: list, sam_paths_disome: list, name_experiment: str,
                 riboseq_assign_at: Union[int, str, list],
                 riboseq_assign_to: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict,
                 footprint_len_experiment: Union[int, list] = None, footprint_len_background: Union[int, list] = None,
                 exclude_genes: list = None, n_core: int = multiprocessing.cpu_count() - 2,
                 verbose=True, recalculate=False):
        super().__init__(temp_repo_dir=temp_repo_dir,
                         sam_paths_background=sam_paths_monosome,
                         sam_paths_experiment=sam_paths_disome,
                         name_experiment=name_experiment,
                         riboseq_assign_at=riboseq_assign_at,
                         riboseq_assign_to=riboseq_assign_to,
                         protein_genome_instance=protein_genome_instance,
                         gene_info_dictionary=gene_info_dictionary,
                         footprint_len_experiment=footprint_len_experiment,
                         footprint_len_background=footprint_len_background,
                         exclude_genes=exclude_genes,
                         n_core=n_core,
                         verbose=verbose,
                         recalculate=recalculate)

        self.output_file_name_fitting_calc = os.path.join(self.temp_repo_dir, f"riboseq_{self.name_experiment}_fitting_calculations.joblib")

        try:
            assert os.access(self.output_file_name_fitting_calc, os.R_OK) and os.path.isfile(self.output_file_name_fitting_calc)
            if self.recalculate:
                print(
                    f"Saved file is found at the path but 'recalculate=True': {self.output_file_name_fitting_calc}.")
                raise AssertionError
            loaded_content = self.load_joblib(self.output_file_name_fitting_calc, self.name_experiment, self.verbose)
            consistent_with = all([
                all([os.path.basename(s) == os.path.basename(l) for s, l in
                     zip(sam_paths_disome, loaded_content["sam_paths_disome"])]),
                all([os.path.basename(s) == os.path.basename(l) for s, l in
                     zip(sam_paths_monosome, loaded_content["sam_paths_monosome"])]),
                os.path.basename(self.output_file_name_fitting_calc) == os.path.basename(loaded_content["output_file_name_fitting_calc"]),
                self.riboseq_assign_at == loaded_content["riboseq_assign_at"],
                len(self.gene_list) == len(loaded_content["gene_list"]),
                all([g == l for g, l in zip(self.gene_list, loaded_content["gene_list"])]),
                # excluded olanlar varsa nolcak
            ])

            if not consistent_with:
                print(
                    f"There is at least one inconsistency between input parameters and loaded content. Recalculating...")
                raise AssertionError
            self.best_model = loaded_content["best_model"]

        except (AssertionError, AttributeError, FileNotFoundError):
            print(f"Sigmoid fitting are being calculated: {self.name_experiment}")
            self.best_model = self.binomial_fitting(self.gene_list, n_core=self.n_core - 2)
            self.save_joblib()

    def save_joblib(self):
        # Write down the output dictionary and list as Joblib object for convenience in later uses.
        if self.verbose:
            print(f"Calculations is being written to directory: {self.temp_repo_dir}")
        to_dump = {"sam_paths_disome": self.sam_paths_experiment,
                   "sam_paths_monosome": self.sam_paths_background,
                   "output_file_name_fitting_calc": self.output_file_name_fitting_calc,
                   "riboseq_assign_at": self.riboseq_assign_at,
                   "gene_list": self.gene_list,
                   "best_model": self.best_model}
        joblib.dump(to_dump, self.output_file_name_fitting_calc)
        if self.verbose:
            print(f"Done: {self.output_file_name_fitting_calc}")

    @staticmethod
    def load_joblib(output_file_name, name_experiment, verbose):
        if verbose:
            print(f"Fitting calculations found for {name_experiment} in path: {output_file_name}")
        return joblib.load(output_file_name)

    def binomial_fitting(self, gene_list, n_core):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if len(gene_list) > 3 * n_core:
                chunk_size = math.ceil(len(gene_list) / n_core)
            else:
                chunk_size = len(gene_list)
            gene_list_chunks = [gene_list[i:i + chunk_size] for i in range(0, len(gene_list), chunk_size)]
            executor = multiprocessing.Pool(len(gene_list_chunks))  # when len(gene_list) > 3 * n_core, len(gene_list_chunks) <= n_core
            result = executor.map(self.binomial_fitting_single_core, gene_list_chunks)
            executor.terminate()
            executor.join()
        result = [j for i in result for j in i]  # Flatten
        return dict(zip(gene_list, result))

    def binomial_fitting_single_core(self, gene_list):
        best_model = list()
        for ind, gene_id in enumerate(gene_list):
            progress_bar(ind, len(gene_list) - 1, verbose=self.verbose)
            best_model.append(self.binomial_fitting_gene(gene_id))
        return best_model

    def create_bf(self, gene_id):  # BF = binomial fitting instance
        # below lines are identical to plot results
        disome_counts = np.sum(self.experiment.gene_assignments[gene_id], axis=0)  # mean problem çıkartıyor,
        monosome_counts = np.sum(self.background.gene_assignments[gene_id], axis=0)
        x_data = np.arange(1, len(disome_counts) + 1)
        return BinomialFitting(x_data, disome_counts, monosome_counts)

    def binomial_fitting_gene(self, gene_id):
        fitter_instance = self.create_bf(gene_id)
        winner_model, model_names, bic_scores, raw_fitting = fitter_instance.select_correct_model()
        return {"winner_model": winner_model,
                "gene_length": self.experiment.gene_assignments[gene_id].shape[1],
                "bic_scores": dict(zip(model_names, bic_scores)),
                "raw_fitting": dict(zip(model_names, raw_fitting))}

    def calculate_curve(self, gene_id, model_name=None):
        results = self.best_model[gene_id]
        if not model_name:
            model_name = results["winner_model"]
        x_data = np.arange(1, results["gene_length"] + 1)
        if model_name == "base":
            return BinomialFitting.model_base(x_data, *results["raw_fitting"]["base"].x)
        elif model_name == "ssig":
            return BinomialFitting.model_ssig(x_data, *results["raw_fitting"]["ssig"].x)
        elif model_name == "dsig":
            return BinomialFitting.model_dsig(x_data, *results["raw_fitting"]["dsig"].x)

    def calculate_onset(self, gene_id, model_name=None):

        def helper(p_ind, p_i, p_j, p_y_mid):
            return p_ind if abs(p_i - p_y_mid) < abs(p_j - p_y_mid) else p_ind + 1

        def find_onset(the_arr, arr_max):
            if len(the_arr) == 0:
                return np.nan
            arr_stride = np.lib.stride_tricks.as_strided(the_arr, shape=(len(the_arr), 2), strides=(the_arr.strides + the_arr.strides))
            # last element is just bullshit
            y_mid = (arr_max + the_arr[0]) / 2
            x_mid = list(set([helper(ind, i, j, y_mid) for ind, (i, j) in enumerate(arr_stride[:-1]) if i <= y_mid <= j]))
            assert len(x_mid) != 0, f"Type 0 error in find_onset: {gene_id}"
            assert len(x_mid) == 1, f"Type 1 error in find_onset: {gene_id}"
            return x_mid[0]

        def find_max(the_arr):
            the_arr = list(the_arr)
            while the_arr:
                max_point = max(the_arr)
                last_point = the_arr.pop(-1)
                if max_point != last_point:
                    return the_arr.index(max_point)
            return 0

        results = self.best_model[gene_id]
        if not model_name:
            model_name = results["winner_model"]

        if model_name == "base":
            return np.nan
        else:
            arr = self.calculate_curve(gene_id, model_name=model_name)
            marr = smooth_array(np.round(arr, 2), window_len=15, window="hanning")  # smooth and round
            derivative = np.gradient(marr)
            descending = np.sum(derivative < 0)
            ascending = np.sum(derivative > 0)
            if model_name == "ssig" and descending != 0:  # negative sigmoid!, only descents
                assert ascending == 0
                return np.nan
            elif model_name == "ssig":
                max_value = arr.max()
                return find_onset(arr, max_value)
            elif model_name == "dsig" and descending != 0 and ascending == 0:  # negative sigmoid!, only descents
                return np.nan
            elif model_name == "dsig" and descending == 0 and ascending != 0:
                max_value = arr.max()
                return find_onset(arr, max_value)
            elif model_name == "dsig":
                max_value_index = find_max(marr)
                sub_arr = arr[:max_value_index]
                max_value = arr[max_value_index]
                return find_onset(sub_arr, max_value)
            else:
                raise AssertionError("Unexpected error!")

    def plot_result(self, gene_id, verbose=False, model_name=None):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            results = self.best_model[gene_id]
            if not model_name:
                model_name = results["winner_model"]
            # below lines are identical to createBD
            f = self.create_bf(gene_id)
            onset = self.calculate_onset(gene_id, model_name=model_name)
            plt.plot(self.calculate_curve(gene_id, model_name="base"), color='black')
            plt.plot(self.calculate_curve(gene_id, model_name="ssig"), color='blue')
            plt.plot(self.calculate_curve(gene_id, model_name="dsig"), color='red')
            plt.vlines(onset, 0, 1, color="blue", alpha=0.75)
            plt.scatter(f.x_data, f.disome / (f.disome + f.monosome), alpha=0.25, color='gray')
            if verbose:
                print(f"base::\nBic: {results['bic_scores']['base']}\nFun: {results['raw_fitting']['base'].fun}\n")
                print(f"ssig::\nBic: {results['bic_scores']['ssig']}\nFun: {results['raw_fitting']['ssig'].fun}\n")
                print(f"dsig::\nBic: {results['bic_scores']['dsig']}\nFun: {results['raw_fitting']['dsig'].fun}\n")
            print(f"Winner: {model_name}\nOnset: {onset}")


class BinomialFitting:

    def __init__(self, x, disome, monosome):
        self.x_data = x  # 1,2,3,4,5...
        self.disome = disome  # 0,0,13,4,0...
        self.monosome = monosome
        assert self.x_data.shape == self.disome.shape == self.monosome.shape, "Shapes do not match. In BinomialFitting."
        self.n_trial = self.disome + self.monosome

    @staticmethod
    def model_base(x, p):
        return p + np.zeros(len(x))

    @staticmethod
    def model_ssig(x, i_init, i_max, a, i_mid):
        return (i_max - i_init) / (1 + np.exp(-a * (x - i_mid))) + i_init

    @staticmethod
    def model_dsig(x, i_init, i_max, i_final, a_1, a_2, i_mid, i_dist):
        return ((i_max - i_init) / (1 + np.exp(-a_1 * (x - i_mid))) + i_init) * \
               ((1 - i_final) / (1 + np.exp(-a_2 * (x - (i_mid + i_dist)))) + i_final)

    # Note: "RuntimeWarning: overflow encountered in exp" will be raised for above two.
    # For most practical purposes, you can probably approximate 1 / (1 + <a large number>) to zero.
    # That is to say, just ignore the warning and move on.
    # Numpy takes care of the approximation for you (when using np.float64).

    # RuntimeWarning: invalid value encountered in subtract df = fun(x) - f0
    # np.nan

    def neg_log_likelihood_base(self, param):
        # The same calculation at model_base(), but for all x_data
        y_predicted = self.model_base(self.x_data, *param)
        negative_log_likelihood = -np.sum(stats.binom.logpmf(k=self.disome, n=self.n_trial, p=y_predicted))
        return negative_log_likelihood

    def neg_log_likelihood_ssig(self, param):
        # The same calculation at model_ssig(), but for all x_data
        y_predicted = self.model_ssig(self.x_data, *param)
        negative_log_likelihood = -np.sum(stats.binom.logpmf(k=self.disome, n=self.n_trial, p=y_predicted))
        return negative_log_likelihood

    def neg_log_likelihood_dsig(self, param):
        # The same calculation at model_dsig(), but for all x_data
        y_predicted = self.model_dsig(self.x_data, *param)
        negative_log_likelihood = -np.sum(stats.binom.logpmf(k=self.disome, n=self.n_trial, p=y_predicted))
        return negative_log_likelihood

    def minimize_base(self, seed=1):
        x0 = np.array([0.5])
        bounds = ((0, 1),)
        minimizer_kwargs = {"method": "SLSQP", "bounds": bounds}
        return basinhopping(func=self.neg_log_likelihood_base, x0=x0, minimizer_kwargs=minimizer_kwargs,
                            seed=seed, niter=150, stepsize=1, interval=50, T=10000)

    def minimize_ssig(self, seed=1):
        x0 = np.array([0.25, 0.65, 0.25, int(len(self.x_data)/2)])
        bounds = ((0, 1), (0, 1), (0, 0.5), (1, len(self.x_data)))
        minimizer_kwargs = {"method": "SLSQP", "bounds": bounds}
        return basinhopping(func=self.neg_log_likelihood_ssig, x0=x0, minimizer_kwargs=minimizer_kwargs,
                            seed=seed, niter=900, stepsize=1, interval=10, T=10000)

    def minimize_dsig(self, seed=1):
        x0 = np.array([0.25, 0.75, 0.5, 0.15, -0.3, int(len(self.x_data) / 2), int(len(self.x_data) / 2)])
        bounds = ((0, 1), (0, 1), (0, 1), (0, 0.5), (-0.5, 0), (1, len(self.x_data)), (1, len(self.x_data)))
        minimizer_kwargs = {"method": "SLSQP", "bounds": bounds}
        return basinhopping(func=self.neg_log_likelihood_dsig, x0=x0, minimizer_kwargs=minimizer_kwargs,
                            seed=seed, niter=900, stepsize=1, interval=10, T=10000)

    def calculate_bic(self, minimized_result):  # Bayesian Information Criterion
        # BIC = -2 * LL + log(N) * k
        # Where log() has the base-e called the natural logarithm, LL is the log-likelihood of the model, N is the number of examples in the training dataset, and k is the number of parameters in the model.
        # The score as defined above is minimized, e.g. the model with the lowest BIC is selected.
        k = len(minimized_result.x)
        ll = -minimized_result.fun
        n = len(self.x_data)
        return -2 * ll + np.log(n) * k

    def select_correct_model(self, seed=1):
        models = ["base", "ssig", "dsig"]
        minimized_results = (self.minimize_base(seed), self.minimize_ssig(seed), self.minimize_dsig(seed))
        bic = [self.calculate_bic(i) for i in minimized_results]
        winner = sorted(zip(bic, models))[0]
        return winner[1], models, bic, minimized_results

    # p = BinomialFitting(np.arange(1,100), np.arange(99), np.ones(99)*100)


class UniprotAnnotation:

    def __init__(self, temp_repo_dir):
        self.download_databases(temp_repo_dir)

    @staticmethod
    def download_databases(temp_repo_dir):
        pass


class ConservationGerp:
    pass


def biomart_mapping(temp_repo_dir, rscript, release, organism):  # organism name should be formatted identical to Organism class
    base_name = os.path.split(os.path.splitext(rscript)[0])[1]
    data_path = os.path.join(temp_repo_dir, f"{base_name}.txt")
    dataset_map = {"homo_sapiens": "hsapiens_gene_ensembl",
                   "mus_musculus": "mmusculus_gene_ensembl"}
    if not os.access(data_path, os.R_OK) or not os.path.isfile(data_path):
        print(f"BiomaRt script is being run for {organism}: {base_name}.")
        r_installation = "RScript" if sys.platform == "darwin" else "Rscript"
        spr = subprocess.run(f"cd {temp_repo_dir}; {which(r_installation)} {rscript} "
                             f"{release} {dataset_map[organism]} {data_path}", shell=True)
        assert spr.returncode == 0, f"Error: {rscript}"

    return pd.read_table(data_path)


# todo: sync this class with pipeline for your conveninence:
class OrganismDatabase:
    ID_MAP = {"homo_sapiens": 9606, "mus_musculus": 10090}

    def __init__(self, organism, ensembl_release, temp_repo_dir):
        assert organism in ["homo_sapiens", "mus_musculus"]
        # todo: also add "homo_sapiens_refseq"

        self.ensembl_release = ensembl_release
        self.organism = organism
        self.organism_id = OrganismDatabase.ID_MAP[self.organism]
        self.temp_repo_dir = temp_repo_dir

        base_temp = f"ftp://ftp.ensembl.org/pub/release-{ensembl_release}"

        if self.organism == "homo_sapiens":
            # Genome GTF
            gtf_temp = f"gtf/homo_sapiens/Homo_sapiens.GRCh38.{self.ensembl_release}.chr_patch_hapl_scaff.gtf.gz"
            self.gtf = os.path.join(base_temp, gtf_temp)
            # Genome GFF3
            gff3_temp = f"gff3/homo_sapiens/Homo_sapiens.GRCh38.{self.ensembl_release}.chr_patch_hapl_scaff.gff3.gz"
            self.gff3 = os.path.join(base_temp, gff3_temp)
            # Genome DNA fasta
            dna_temp = "fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
            self.dna = os.path.join(base_temp, dna_temp)
            # Transcriptome DNA fasta
            cdna_temp = "fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
            self.cdna = os.path.join(base_temp, cdna_temp)
            # Protein Fasta
            pep_temp = "fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
            self.pep = os.path.join(base_temp, pep_temp)

        elif self.organism == "mus_musculus":
            # Genome GTF
            gtf_temp = f"gtf/mus_musculus/Mus_musculus.GRCm38.{self.ensembl_release}.chr_patch_hapl_scaff.gtf.gz"
            self.gtf = os.path.join(base_temp, gtf_temp)
            # Genome GFF3
            gff3_temp = f"gff3/mus_musculus/Mus_musculus.GRCm38.{self.ensembl_release}.chr_patch_hapl_scaff.gff3.gz"
            self.gff3 = os.path.join(base_temp, gff3_temp)
            # Genome DNA fasta
            dna_temp = "fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
            self.dna = os.path.join(base_temp, dna_temp)
            # Transcriptome DNA fasta
            cdna_temp = "fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
            self.cdna = os.path.join(base_temp, cdna_temp)
            # Protein Fasta
            pep_temp = "fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz"
            self.pep = os.path.join(base_temp, pep_temp)

    def get_db(self, db):
        db_url = eval(f"self.{db}")
        output_path_compressed = os.path.join(self.temp_repo_dir, os.path.basename(db_url))
        output_path_uncompressed = os.path.splitext(output_path_compressed)[0]
        if not os.access(output_path_uncompressed, os.R_OK) or not os.path.isfile(output_path_uncompressed):
            print(f"Downloading from the server for {db}:{os.linesep}{db_url}")
            if not os.access(output_path_compressed, os.R_OK) or not os.path.isfile(output_path_compressed):
                subprocess.run(f"cd {self.temp_repo_dir}; curl -L -O --silent {db_url}", shell=True)
            subprocess.run(f"cd {self.temp_repo_dir}; gzip -d -q {output_path_compressed}", shell=True)
        return output_path_uncompressed


def ensembl_release_object_creator(temp_repo_dir, release, organism):

    organism_instance = OrganismDatabase(organism, release, temp_repo_dir)
    gtf_path = organism_instance.get_db("gtf")
    transcript_fasta_path = organism_instance.get_db("cdna")
    protein_fasta_path = organism_instance.get_db("pep")

    ero = pyensembl.Genome(
        reference_name=f"{organism}_{release}",
        annotation_version=release,
        annotation_name=f"{organism}_{release}",
        gtf_path_or_url=gtf_path,
        transcript_fasta_paths_or_urls=[transcript_fasta_path],
        protein_fasta_paths_or_urls=[protein_fasta_path],
        copy_local_files_to_cache=False,
        cache_directory_path=temp_repo_dir
    )
    ero.index()
    return ero


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


def gene_entire_cds(protein_genome_instance, gene_info_dictionary, gene_id, make_array=False):
    transcript_list = gene_info_dictionary[gene_id].transcripts["ensembl_transcript_id"].to_list()
    cds_all_ranges = [protein_genome_instance.db[transcript_id][0] for transcript_id in transcript_list]
    cds_all_ranges = reduce_range_list([j for i in cds_all_ranges for j in i])  # it sorts
    strand = -1 if protein_genome_instance.db[transcript_list[0]][4] == "-" else 1  # all tx are always the same
    # Ters olanları ters sırada veriyor!!!
    # ProteinGenome zaten düzeltiyor transcript'i eğer -1'deyse ama reduce_range_list sırayı bozuyor.
    if strand == -1:
        cds_all_ranges = [[j, i] for i, j in reversed(cds_all_ranges)]  # reverse if at '-'
    if not make_array:
        return cds_all_ranges
    else:
        return np.concatenate([np.arange(i, j + strand, strand) for i, j in cds_all_ranges])


def best_transcript_cds(protein_genome_instance, gene_info_dictionary, gene_id, make_array=False):
    best_transcript = gene_info_dictionary[gene_id].transcripts.iloc[0][0]  # At least 1 transcript exists
    cds_ranges = protein_genome_instance.db[best_transcript][0]
    # Ters olanları ters sırada veriyor!!!
    # ProteinGenome zaten düzeltiyor transcript'i eğer -1'deyse
    if not make_array:
        return cds_ranges
    else:
        strand = -1 if protein_genome_instance.db[best_transcript][4] == "-" else 1
        return np.concatenate([np.arange(i, j + strand, strand) for i, j in cds_ranges])


def ensembl_range_sum(ranges: list) -> int:
    """
    The function is to calculate the total CDS/exon length of a transcript.
    :param ranges: Nested list. Obtained by coding_sequence_position_ranges or exon_intervals method of tx object
    :return: Total length of the CDS/exon
    """
    total = 0  # Initialize
    for i, j in ranges:  # For each CDS/exon of a transcript, get start and end positions
        # j always larger or equal to i
        total += (j - i + 1)  # Plus 1 is due to the fact that start and end is included.
    return total


def array_to_ranges(i): # set or list
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


def reduce_range_list(ranges):
    rgs = sorted(ranges, key=lambda pair: (pair[0], - pair[1]))
    arr = set([k for i, j in rgs for k in range(i, j + 1)])
    return list(array_to_ranges(arr))


def save_joblib(object_to_dump, absolute_output_path, verbose):
    # Write down the output dictionary and list as Joblib object for convenience in later uses.
    print_verbose(verbose, f"Instance is being written to directory.")
    joblib.dump(object_to_dump, absolute_output_path)
    print_verbose(verbose, f"Done: {absolute_output_path}")


def load_joblib(object_name, output_file_name, verbose):
    print_verbose(verbose, f"{object_name} found in path: {output_file_name}")
    return joblib.load(output_file_name)


def progress_bar(iteration: int, total: int, prefix: str = 'Progress:', suffix: str = '', decimals: int = 1,
                 bar_length: int = 20, verbose: bool = False):
    """
    This function should be called inside of loop, gives the loop's progress.
    :param iteration: Current iteration.
    :param total: Total iteration.
    :param prefix: Placed before progress bar.
    :param suffix: Placed after progress bar.
    :param decimals: Number of decimals in percent complete.
    :param bar_length: Character length of bar.
    :param verbose: For convenience in some uses
    :return: Void function. Nothing is returned.
    """
    if verbose:
        filled_length = int(round(bar_length * iteration / float(total)))
        percents = round(100.00 * (iteration / float(total)), decimals)
        bar = '█' * filled_length + '-' * (bar_length - filled_length)
        sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
        sys.stdout.flush()
        if iteration == total:
            sys.stdout.write('\n')
            sys.stdout.flush()


def smooth_array(x: np.ndarray, window_len: int, window: str) -> np.ndarray:
    """
    Smooth the data using a window with requested size. This method is based on the convolution of a scaled window
    with the signal. The signal is prepared by introducing reflected copies of the signal (with the window size) in
    both ends so that transient parts are minimized in the beginning and end part of the output signal.
    Adapted from: http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    :param x: The input signal
    :param window_len: The dimension of the smoothing window.
    :param window: The type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'.
    Note that flat window will produce a moving average smoothing.
    :return: The smoothed signal
    """
    # Check if the parameters are set correctly.
    assert x.ndim == 1, "Error: smooth() only accepts 1 dimension arrays."
    assert x.size >= window_len, "Error: Input vector needs to be bigger than window size."
    assert window in ('flat', 'hanning', 'hamming', 'bartlett', 'blackman'), "Error: Window type is unknown"
    assert window_len % 2 == 1, "Window length should be odd integer."
    # No need for further calculation if window length is smaller than 3
    if window_len < 3:
        return x
    # Create the kernel for convolution and calculate the convolution
    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    w = np.ones(window_len, 'd') if window == 'flat' else eval('np.' + window + '(window_len)')
    y = np.convolve(w / w.sum(), s, mode='valid')
    # Correct the shape of the result
    ds = y.shape[0] - x.shape[0]  # difference of shape
    dsb = ds // 2  # [almost] half of the difference of shape for indexing at the beginning
    dse = ds - dsb  # rest of the difference of shape for indexing at the end
    y = y[dsb:-dse]
    # Return the result
    return y


def print_verbose(verbose, message):
    if verbose:
        print(f"{Col.G}[{time.strftime('%d/%m/%Y %H:%M:%S %Z')}]{Col.E} {message}")


def cleartarget(dir_path: str):
    """
    Delete everything in the output folder. Creates new one with the same name
    :param dir_path: Path of directory to delete
    """
    try:
        rmtree(dir_path)
    except FileNotFoundError:
        pass
    os.mkdir(dir_path)


class Col:
    H = '\033[95m\033[1m'  # Header
    B = '\033[94m'  # Blue
    C = '\033[96m'  # Cyan
    G = '\033[92m'  # Green
    W = '\033[93m'  # Warning
    F = '\033[91m'  # Fail
    E = '\033[0m'  # End
    D = '\033[1m'  # Bold
    U = '\033[4m'  # Underline


# ilia said: t1>t2 always
# ilia said: instead of np.sum calculate twice, get onsets, get average


# End
