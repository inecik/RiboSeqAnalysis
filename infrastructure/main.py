"""

"""


import sys
import os
import re
import subprocess
import sys
import warnings
import itertools
# import logging  # todo

import joblib
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing
from scipy.optimize import basinhopping
# from statsmodels.stats.proportion import proportion_confint  # todo
import pyensembl
import math
import pysam
from shutil import which
from Bio.Seq import translate
from Bio import pairwise2
import matplotlib.pyplot as plt
import random


script_path_infrastructure = __file__


class Infrastructre:

    def __init__(self, temp_repo_dir, exclude_genes=list(), ensembl_release=102,
                 sixtymers=None, serb=None, coco=None, riboseq_assign_to="best_transcript", riboseq_assign_at=-1,
                 include_gene3d=False, verbose=True):

        self.temp_repo_dir = temp_repo_dir
        self.exclude_genes = exclude_genes
        self.ensembl_release = ensembl_release
        self.verbose = verbose

        self.script_path_infrastructure = script_path_infrastructure
        self.script_path_gene_info_db = os.path.abspath(os.path.join(os.path.dirname(self.script_path_infrastructure), "gene_info_database.R"))
        self.script_path_gene_info_names_db = os.path.abspath(os.path.join(os.path.dirname(self.script_path_infrastructure), "gene_info_names.R"))
        self.script_path_gene_info_uniprot_db = os.path.abspath(os.path.join(os.path.dirname(self.script_path_infrastructure), "gene_info_uniprot.R"))

        gene_info_database = biomart_mapping(self.temp_repo_dir, self.script_path_gene_info_db, self.ensembl_release)
        gene_info_names = biomart_mapping(self.temp_repo_dir, self.script_path_gene_info_names_db, self.ensembl_release)
        gene_info_uniprot = biomart_mapping(self.temp_repo_dir, self.script_path_gene_info_uniprot_db, self.ensembl_release)

        gene_info_database = gene_info_database[gene_info_database["transcript_biotype"] == "protein_coding"]
        self.gene_list = sorted(np.unique(gene_info_database["ensembl_gene_id"].dropna()))
        self.gene_info = gene_class_dict_generate(self.temp_repo_dir, self.gene_list, gene_info_database, gene_info_names, gene_info_uniprot, verbose=self.verbose)

        ero = ensembl_release_object_creator(self.temp_repo_dir, self.ensembl_release)

        transcript_list = sorted(np.unique(gene_info_database["ensembl_transcript_id"].dropna()))
        self.protein_genome = ProteinGenome(self.temp_repo_dir, transcript_list, ero, verbose=self.verbose)

        # Integrate RiboSeq data

        self.riboseq_assign_to = riboseq_assign_to
        self.riboseq_assign_at = riboseq_assign_at

        # Should be first to calculate, since multiprocessing is quite memory inefficient.
        if coco:  # [monosome_sam, disome_sam]
            self.riboseq_coco = RiboSeqCoco(self.temp_repo_dir, coco[0], coco[1], self.riboseq_assign_at,
                                            self.riboseq_assign_to, self.protein_genome, self.gene_info,
                                            exclude_genes=self.exclude_genes, verbose=self.verbose)
        else:
            self.riboseq_coco = None

        if sixtymers:  # [translatome_sam, sixtymer_sam]
            self.riboseq_sixtymers = RiboSeqSixtymers(self.temp_repo_dir, sixtymers[0], sixtymers[1], self.riboseq_assign_at,
                                                      self.riboseq_assign_to, self.protein_genome, self.gene_info, exclude_genes=self.exclude_genes, verbose=self.verbose)
        else:
            self.riboseq_sixtymers = None

        if serb:  # [[str:name, str:translatome_sam, str:experiment_sam], [...]]
            self.riboseq_serb = dict()
            for serb0, serb1, serb2 in serb:
                self.riboseq_serb[serb0] = RiboSeqSelective(self.temp_repo_dir, serb1, serb2, serb0, self.riboseq_assign_at,
                                                             self.riboseq_assign_to, self.protein_genome, self.gene_info, exclude_genes=self.exclude_genes, verbose=self.verbose)
            else:
                self.riboseq_serb = None

        # Integrate protein annotations

        if include_gene3d:
            self.script_path_gene3d_db = os.path.abspath(os.path.join(os.path.dirname(self.script_path_infrastructure), "gene3d.R"))
            self.gene3d_database = EnsemblDomain(self.temp_repo_dir, self.script_path_gene3d_db, self.protein_genome, ero, verbose=self.verbose)



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
    :param verbose:  If True, it will print to stdout about the process computer currently calculates.
    :return: Dictionary of gene_id → Gene instance
    """

    # Get the absolute output file name.
    output_path = os.path.join(temp_repo_dir, "gene_info_database.joblib")
    # Check if there is already a calculated object saved before.
    if not os.access(output_path, os.R_OK) or not os.path.isfile(output_path) or recalculate:
        # If no file exists, or overwrite is 'True'.
        output = dict()  # Initiate an dictionary to fill up.
        print_verbose(verbose, f"{Col.HEADER}Gene information dictionary are being created.{Col.ENDC}")
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
        print_verbose(verbose, f"{Col.HEADER}Results are being written to directory: {temp_repo_dir}{Col.ENDC}")
        # Save the resulting filled dictionary into a joblib file to load without calculating again in next runs.
        joblib.dump(output, output_path)
        print_verbose(verbose, f"Done: {output_path}")
        return output  # Return the resulting dictionary.
    else:  # If the joblib file is found in the temp directory, or overwrite is 'False'.
        print_verbose(verbose, f"{Col.HEADER}Gene information dictionary is found in path: {output_path}{Col.ENDC}")
        return joblib.load(output_path)  # Load and return the resulting dictionary.


class Gene:
    """
    Stores all the relevant data for a given gene, prioritizes the transcripts.
    """

    def __init__(self, one_gene_df: pd.DataFrame, one_gene_names: pd.DataFrame, one_gene_uniprot: pd.DataFrame):
        """
        'A=k[k["ensembl_gene_id"] == gene_id]' where k is biomart_mapping() output.
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
        assert len(np.unique(transcripts["transcript_mane_select"].dropna())) <= 1
        # Based on my research, sorting the transcripts by the following order will always give the best possible
        # transcript in terms of functional relevance and amount of evidence supporting its presence.
        transcripts.sort_values(by=["transcript_mane_select", "transcript_appris", "transcript_gencode_basic",
                                    "transcript_tsl", "external_transcript_name"], inplace=True, ignore_index=True)
        # See https://m.ensembl.org/info/genome/genebuild/transcript_quality_tags.html for more information
        return transcripts  # Return the sorted and narrowed dataframe


class ProteinGenome:
    """
    Store all transcripts sequence as well as corresponding proteins' sequence. This class is created to solve out some
    problems in the databases that cause my calculations to suffer. For example, for some transcripts that has
    undefined 3' or 5', the length of protein sequence does not match the length of transcript sequence. The methods in
    this class solves the problem by adding extra Ns to the transcript sequence. Note that there is no example (in
    Ensembl version 102) that has a transcript length which is longer than protein length (when multiplied with 3, of
    course); the issue is in the opposite way. In addition to this, the class convert genomic coordinates to protein or
    transcript ranges and vice versa.
    """

    def __init__(self, temp_repo_dir: str, transcript_list: list, ensembl_release_object: pyensembl.Genome,
                 verbose: bool = True, recalculate: bool = False):
        """
        :param temp_repo_dir: Full path directory where temporary files will be stored.
        :param transcript_list: List of gene IDs, for which to create the dictionary. Ideally, it should be sorted.
        :param ensembl_release_object: Output of ensembl_release_object_creator()
        :param verbose:  If True, it will print to stdout about the process computer currently calculates.
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
                print_verbose(verbose, f"{Col.WARNING}Saved file is found at the path but 'recalculate' is activated: "
                                       f"{self.output_file_name}{Col.ENDC}.")
                raise AssertionError  # Raise the error to go to the except statement.
            # Load the saved content from the directory.
            loaded_content = self.load_joblib(self.output_file_name, self.verbose)
            # Check the current run is consistent with the saved file
            consistent_with = all([  # Make sure all below is True
                self.transcript_list == loaded_content.transcript_list,  # All transcripts are the same, and sorted.
                self.ensembl_release == loaded_content.ensembl_release,  # The same ensembl release is used.
            ])
            if not consistent_with:  # If there is a problem, the stored does not match with current run.
                print_verbose(verbose, f"{Col.WARNING}There is at least one inconsistency between input parameters and "
                                       f"loaded content. Recalculating.{Col.ENDC}")
                raise AssertionError  # Raise the error to go to the except statement.
            # Otherwise, just accept the database saved in previous run.
            self.db = loaded_content.db
        except (AssertionError, AttributeError, FileNotFoundError):  # If an error is raised.
            print_verbose(verbose, f"{Col.HEADER}Protein genome mapping are being calculated.{Col.ENDC}")
            self.db = self.calculate_transcript_mapping(ensembl_release_object)  # Calculate to get the mappings
            # Save the resulting filled dictionary into a joblib file to load without calculating again in next runs.
            self.save_joblib()

    @staticmethod
    def consistent_coding_ranges(transcript_object: pyensembl.Transcript):
        """
        00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        :param transcript_object:
        :return:
        """

        def protein_coding_ranges(transcript_object):
            """
            Some transcript has no defined start or stop codon. This causes the pyensembl module to fail in to find
            coding sequence with coding_sequence method. This function is to replace this method.
            Tested: for a given transcript, all CDSs are at the same strand.
            Tested: for a given transcript, coding_sequence_position_ranges method gives CDSs in order.
            Tested: for transcripts which does not raise error with coding_sequence, this function gives identical results.
            :param transcript_object: Created by ensembl_release_object_creator() function
            :return: List for Coding sequence. Sequence: "".join([transcript_object.sequence[s: e + 1] for s, e in cdr_relative])
            """
            transcript_object_exons = transcript_object.exons
            exon_numbers = transcript_object.db.query(select_column_names=["exon_number"], feature="CDS",
                                                      filter_column="transcript_id", filter_value=transcript_object.id)
            first_coding_exon = int(exon_numbers[0][0]) - 1
            origin_exon = transcript_object_exons[first_coding_exon]
            offset = ensembl_range_sum(transcript_object.exon_intervals[:first_coding_exon])

            cdr = transcript_object.coding_sequence_position_ranges
            if transcript_object.strand == "+":
                cdr_relative_temp = [[s - origin_exon.start, e - origin_exon.start] for s, e in cdr]
            else:  # Notice: it also switches start and end positions of the ranges
                cdr_relative_temp = [[origin_exon.end - e, origin_exon.end - s] for s, e in cdr]

            cdr_relative = [cdr_relative_temp.pop(0)]
            while len(cdr_relative_temp) > 0:
                intron = cdr_relative_temp[0][0] - cdr_relative[-1][1] - 1
                cdr_relative_temp = [[s - intron, e - intron] for s, e in cdr_relative_temp]
                cdr_relative.append(cdr_relative_temp.pop(0))

            cdr_relative = [[offset + s, offset + e] for s, e in cdr_relative]  # Bring back skipped exons

            return cdr_relative

        def transcript_frame(query_transcript, target_protein, table):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                results = list()
                for frame_start in range(3):
                    frame = translate(query_transcript[frame_start:], table=table)
                    score = pairwise2.align.localxx(frame, target_protein, score_only=True)
                    if score:
                        results.append([score, frame_start])
                assert len(results) > 0, "No alignment"
                return sorted(results, reverse=True)[0][1]

        def find_overhang(query, target, reverse, search_in=10):  # arbitrary 10
            half_search = math.floor(search_in / 2)
            target_temp = target[-search_in:][::-1] if reverse else target[:search_in]
            query_temp = query[-search_in:][::-1] if reverse else query[:search_in]
            for query_offset in list(range(0, half_search)):
                query_last = query_temp[query_offset: query_offset + half_search]
                relative_position = target_temp.find(query_last)
                if relative_position != -1:
                    return relative_position - query_offset
            raise AssertionError("Error in find_overhang")

        try:
            coding_position_ranges = [list(i) for i in transcript_object.coding_sequence_position_ranges]
        except ValueError:  # ValueError: Transcript does not contain feature CDS
            return [[i, j] if transcript_object.strand == "+" else [j, i] for i, j in transcript_object.exon_intervals], \
                   None, transcript_object.sequence, transcript_object.contig, transcript_object.strand

        if transcript_object.complete:
            assert len(transcript_object.coding_sequence) == len(transcript_object.protein_sequence) * 3 + 3
            assert ensembl_range_sum(coding_position_ranges) == len(transcript_object.protein_sequence) * 3
            return [[i, j] if transcript_object.strand == "+" else [j, i] for i, j in coding_position_ranges], \
                   transcript_object.protein_sequence, transcript_object.coding_sequence[:-3], transcript_object.contig, transcript_object.strand
        else:
            target = transcript_object.protein_sequence
            query_positions = protein_coding_ranges(transcript_object)
            query_transcript = "".join([transcript_object.sequence[s: e + 1] for s, e in query_positions])
            translate_table = 2 if transcript_object.contig == "MT" else 1
            frame_start = transcript_frame(query_transcript, target, translate_table)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                query = translate(query_transcript[frame_start:], table=translate_table)
            if query[-1] == "*":
                query = query[:-1]  # Selanosistein falan * oluyor görünüyor çünkü
            five_overhang = find_overhang(query, target, reverse=False)
            three_overhang = find_overhang(query, target, reverse=True)
            three_overhang_additional = ensembl_range_sum(coding_position_ranges) - frame_start - 3 * len(query)

            mofify5 = -3 * five_overhang + frame_start
            modify3 = 3 * three_overhang - three_overhang_additional

            if transcript_object.strand == "+":
                coding_position_ranges[0][0] += mofify5
                coding_position_ranges[-1][1] += modify3
            else:
                coding_position_ranges[0][1] -= mofify5
                coding_position_ranges[-1][0] -= modify3

            query_transcript = query_transcript[-mofify5:] if mofify5 >= 0 else "N" * -mofify5 + query_transcript
            query_transcript = query_transcript + "N" * modify3 if modify3 >= 0 else query_transcript[:modify3]

            assert ensembl_range_sum(coding_position_ranges) == len(transcript_object.protein_sequence) * 3
            assert len(query_transcript) == len(transcript_object.protein_sequence) * 3

            return [[i, j] if transcript_object.strand == "+" else [j, i] for i, j in coding_position_ranges], \
                   transcript_object.protein_sequence, query_transcript, transcript_object.contig, transcript_object.strand

    def calculate_transcript_mapping(self, ensembl_release_object):
        assert len(np.unique(self.transcript_list)) == len(self.transcript_list)
        output = dict()
        for ind, transcript_id in enumerate(self.transcript_list):
            if self.verbose:
                progress_bar(ind, len(self.transcript_list) - 1)
            transcript_object = ensembl_release_object.transcript_by_id(transcript_id)
            output[transcript_id] = self.consistent_coding_ranges(transcript_object)
        return output

    def save_joblib(self):
        # Write down the output dictionary and list as Joblib object for convenience in later uses.
        if self.verbose:
            print(f"{Col.HEADER}Instance is being written to directory: {self.temp_repo_dir}{Col.ENDC}")
        joblib.dump(self, self.output_file_name)
        if self.verbose:
            print(f"Done: {self.output_file_name}")

    @staticmethod
    def load_joblib(output_file_name, verbose):
        if verbose:
            print(f"{Col.HEADER}ProteinGenome instance found in path: {output_file_name}{Col.ENDC}")
        return joblib.load(output_file_name)

    def protein2genome(self, transcript_id, start, end):
        transcript = self.db[transcript_id]

        assert 1 <= start <= end <= len(transcript[1]), f"Wrong range for transcript {transcript_id}: min: 1, max: {len(transcript[1])}"
        start_nt = start * 3 - 3  # -1: convert python index
        end_nt = end * 3 - 1

        origin = transcript[0][0][0]
        cdr_relative_temp = [[i - origin, j - origin] if transcript[4] == "+" else [origin - i, origin - j] for i, j in transcript[0]]
        cdr_relative = [cdr_relative_temp.pop(0)]
        introns = [0]
        while len(cdr_relative_temp) > 0:
            intron = cdr_relative_temp[0][0] - cdr_relative[-1][1] - 1
            cdr_relative_temp = [[s - intron, e - intron] for s, e in cdr_relative_temp]
            cdr_relative.append(cdr_relative_temp.pop(0))
            introns.append(intron)

        introns = list(itertools.accumulate(introns))
        output_relative_introns = list()
        for ind, (cdr_start, cdr_end) in enumerate(cdr_relative):
            if start_nt > cdr_end:
                continue
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

        return [[i + origin, j + origin] if transcript[4] == "+" else [origin - i, origin - j] for i, j in output_relative_introns]

    @staticmethod
    def genome_to_sequences(self, chromosome, ranges):
        pass


class EnsemblDomain:

    def __init__(self, temp_repo_dir, rscript, protein_genome_instance=None, ensembl_release_object=None, ensembl_release=102, verbose=True, recalculate=False):
        self.temp_repo_dir = temp_repo_dir
        self.rscript = rscript
        self.base_name = os.path.split(os.path.splitext(self.rscript)[0])[1]
        self.data_path = os.path.join(self.temp_repo_dir, f"{self.base_name}.txt")
        self.output_file_name = os.path.join(self.temp_repo_dir, f"{self.base_name}.joblib")
        self.ensembl_release = ensembl_release
        self.verbose = verbose
        self.recalculate = recalculate

        try:
            assert os.access(self.output_file_name, os.R_OK) and os.path.isfile(self.output_file_name)
            if self.recalculate:
                print(f"{Col.WARNING}Saved file is found at the path but 'recalculate=True': {self.output_file_name}{Col.ENDC}.")
                raise AssertionError
            loaded_content = self.load_joblib(self.output_file_name, self.base_name, self.verbose)
            consistent_with = all([
                self.base_name == loaded_content.base_name,
                self.ensembl_release == loaded_content.ensembl_release
            ])
            if not consistent_with:
                print(f"{Col.WARNING}There is at least one inconsistency between input parameters and loaded content. Recalculating...{Col.ENDC}")
                raise AssertionError
            self.df = loaded_content.df
            self.columns = loaded_content.columns
        except (AssertionError, AttributeError, FileNotFoundError):
            assert protein_genome_instance, "Undefined protein_genome_instance."
            assert ensembl_release_object, "Undefined ensembl_release_object."
            print(f"{Col.HEADER}Ensembl domains are being calculated: {self.base_name}{Col.ENDC}")
            self.df = biomart_mapping(self.temp_repo_dir, self.rscript, self.ensembl_release)
            self.columns = self.df.columns
            self.convert2genome(protein_genome_instance, ensembl_release_object)
            self.save_joblib()

    def convert2genome(self, protein_genome_instance, ensembl_release_object):
        protein_ids = self.df[self.columns[0]]
        domain_starts = self.df[self.columns[1]]
        domain_ends = self.df[self.columns[2]]
        coordinates_ranges = list()
        coordinates_contig = list()
        for ind, (protein_id, domain_start, domain_end) in enumerate(zip(protein_ids, domain_starts, domain_ends)):
            if self.verbose:
                progress_bar(ind, len(protein_ids) - 1)
            transcript_id = ensembl_release_object.transcript_id_of_protein_id(protein_id)
            # orientation doğru mu???
            coordinates_ranges.append(protein_genome_instance.protein2genome(transcript_id, domain_start, domain_end))
            coordinates_contig.append(protein_genome_instance.db[transcript_id][3])
        self.df["genome_chromosome"] = coordinates_contig
        self.df["genome_coordinate"] = coordinates_ranges

    def save_joblib(self):
        # Write down the output dictionary and list as Joblib object for convenience in later uses.
        if self.verbose:
            print(f"{Col.HEADER}Instance is being written to directory: {self.temp_repo_dir}{Col.ENDC}")
        joblib.dump(self, self.output_file_name)
        if self.verbose:
            print(f"Done: {self.output_file_name}")

    @staticmethod
    def load_joblib(output_file_name, base_name, verbose):
        if verbose:
            print(f"{Col.HEADER}Ensembl database instance found for {base_name} in path: {output_file_name}{Col.ENDC}")
        return joblib.load(output_file_name)

# TODO: co-co site'ı hesaplayan fitting şeyi yaz: coco_site'ı eklemlendir.
#   fitting'de gerçekten raw sayılar mı yoksa rpkm falan mı?
#   makaledeki şeyleri eksiksiz eklemlendir (rpkm, ci, rolling window etc.)
#   lowCI highCI falan hesapla, kolaysa metagene'i de eklemlendir değilse boşver
#   arpat yöntemine göz at, doğru yazıldığından emin ol
#   birkaç peak detection tool'u ekle
#   conservation'u ekle
#   gene class'ını bitir.
#   uniprot'u eklemeden önce infrastructure'ı bitir ve buralara
#   bol bol açıklama ekle. bütün koda (protein_sequence asıl: nedenini açıkla)

# TODO: SORU: Sadece protein coding mi, yoksa her şey mi?
# şu anda bütün genleri alıyor rrna falan dahil, bu rpkm'i yapay olarak arttıracaktır..
# alignment'da değil sadece burada protein-coding only demeliyiz bence. en başta infrastructure'ı oluştururken

class RiboSeqAssignment:

    def __init__(self, sam_paths: list, temp_repo_dir: str, assignment: int, selection: str, riboseq_group: str,
                 protein_genome_instance: ProteinGenome, gene_info_dictionary: dict, verbose=True, recalculate=False):
        # :param riboseq_group: String to identify the RiboSeq experiment annotation; like "translatome" or "60mers"
        self.sam_paths = sorted(sam_paths)
        self.temp_repo_dir = temp_repo_dir
        self.output_file_name = os.path.join(temp_repo_dir, f"riboseq_{riboseq_group}_on_{selection}.joblib")
        self.gene_list = sorted(gene_info_dictionary.keys())
        self.exclude_gene_list = list()
        self.assignment = assignment
        self.verbose = verbose
        self.selection = selection
        self.riboseq_group = riboseq_group  # Define riboseq_group variable for a function.
        self.recalculate = recalculate

        try:
            assert os.access(self.output_file_name, os.R_OK) and os.path.isfile(self.output_file_name)
            if self.recalculate:
                print(f"{Col.WARNING}Saved file is found at the path but 'recalculate=True': {self.output_file_name}{Col.ENDC}.")
                raise AssertionError
            loaded_content = self.load_joblib(self.riboseq_group, self.output_file_name, self.verbose)
            consistent_with = all([
                all([os.path.basename(s) == os.path.basename(l) for s, l in zip(self.sam_paths, loaded_content.sam_paths)]),
                os.path.basename(self.output_file_name) == os.path.basename(loaded_content.output_file_name),
                self.assignment == loaded_content.assignment,
                len(self.gene_list) == len(loaded_content.gene_list),
                all([g == l for g, l in zip(self.gene_list, loaded_content.gene_list)]),
            ])
            if not consistent_with:
                print(f"{Col.WARNING}There is at least one inconsistency between input parameters and loaded content. Recalculating...{Col.ENDC}")
                raise AssertionError
            self.gene_assignments = loaded_content.gene_assignments
            self.gene_lengths = loaded_content.gene_lengths
            self.total_assigned_gene = loaded_content.total_assigned_gene
            self.total_assigned = loaded_content.total_assigned

        except (AssertionError, AttributeError, FileNotFoundError):
            print(f"{Col.HEADER}Gene assignments are being calculated: {self.riboseq_group}{Col.ENDC}")
            self.calculate_gene_assignments(protein_genome_instance, gene_info_dictionary)
            self.gene_lengths = dict(zip(self.gene_list, [int(self.gene_assignments[i].shape[1]) for i in self.gene_list]))
            self.exclude_genes_calculate_stats(self.exclude_gene_list)
            self.save_joblib()

    def exclude_genes_calculate_stats(self, exclude_gene_list: list):
        for gene_id in exclude_gene_list:
            target_shape = self.gene_assignments[gene_id].shape
            empty_matrix = np.empty(target_shape)
            empty_matrix.fill(np.nan)
            self.gene_assignments[gene_id] = empty_matrix

        self.exclude_gene_list.extend(exclude_gene_list)
        self.total_assigned_gene = dict(zip(self.gene_list, [np.nansum(self.gene_assignments[i], axis=1) for i in self.gene_list]))
        self.total_assigned = np.nansum(np.array(list(self.total_assigned_gene.values())), axis=0)

    def calculate_gene_assignments(self, protein_genome_instance, gene_info_dictionary):

        if self.verbose:
            print(f"{Col.HEADER}Genes are allocated to chromosomes.{Col.ENDC}")
        # Get the mapping from chromosomes to list of genes which they carry
        chromosome_gene = self.chromosomes_genes_matcher(gene_info_dictionary, self.gene_list)
        if self.verbose:
            print(f"{Col.HEADER}CDS ranges of genes are being calculated.{Col.ENDC}")
        # Get array of genomic positions for all genes, create a dictionary out of it.
        if self.selection == "gene":
            positions_gene = {gene_id: gene_entire_cds(protein_genome_instance, gene_info_dictionary, gene_id, make_array=True) for gene_id in self.gene_list}
        elif self.selection == "best_transcript":
            positions_gene = {gene_id: best_transcript_cds(protein_genome_instance, gene_info_dictionary, gene_id, make_array=True) for gene_id in self.gene_list}
        else:
            raise AssertionError("Selection variable must be either 'gene' or 'best_transcript'.")

        if self.verbose:
            print(f"{Col.HEADER}Footprint are being assigned to genomic coordinates.{Col.ENDC}")
        # Assign footprints to genomic positions
        footprint_genome_assignment_list = [self.footprint_assignment(sam_path, assignment=self.assignment, verbose=self.verbose) for sam_path in self.sam_paths]

        if self.verbose:
            print(f"{Col.HEADER}Footprint counts are being calculated and assigned to genes.{Col.ENDC}")
        # Assign footprints to gene positions
        self.gene_assignments = self.footprint_counts_to_genes(footprint_genome_assignment_list, chromosome_gene, positions_gene, verbose=self.verbose)

    def save_joblib(self):
        # Write down the output dictionary and list as Joblib object for convenience in later uses.
        if self.verbose:
            print(f"{Col.HEADER}Instance is being written to directory: {self.temp_repo_dir}{Col.ENDC}")
        joblib.dump(self, self.output_file_name)
        if self.verbose:
            print(f"Done: {self.output_file_name}")

    @staticmethod
    def load_joblib(riboseq_group, output_file_name, verbose):
        if verbose:
            print(f"{Col.HEADER}RiboSeq assignment found for {riboseq_group} in path: {output_file_name}{Col.ENDC}")
        return joblib.load(output_file_name)

    @staticmethod
    def chromosomes_genes_matcher(gene_info_dictionary, gene_list):
        """
        This script is to assign genes into chromosomes.
        :param gene_info_dictionary: Pass
        :param gene_list: List of genes to assign chromosomes
        :return: Dictionary of assigned genes.
        """
        chromosomes = dict()  # Initialize a dictionary to fill up
        for gene_id in gene_list:  # For each gene in the gene_list
            contig = gene_info_dictionary[gene_id].chromosome
            # Find the chromosome name (contig), add the gene into the relevant list
            chromosomes[contig] = chromosomes.get(contig, []) + [gene_id]
        return chromosomes

    @staticmethod
    def footprint_assignment(sam_path, assignment, verbose=False):
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

        counter = 0
        with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open SAM file with pysam library
            sam_iterator = sam_handle.fetch()  # Get the iterator to iterate through reads in for loop
            assigned_positions = dict()  # Initialize dictionary with keys composed of chromosome names

            for ind, e in enumerate(sam_iterator):  # Iterate through the entries

                # If 'verbose' activated, print the progress bar in a meaningful intervals
                if verbose and (ind % 1000 == 0 or ind == iteration - 1):
                    progress_bar(ind, iteration - 1)

                assert sum(e.get_cigar_stats()[1][5:]) == 0, "Make sure cigar string composed of D, I, M, N, S only!"
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
                    assert not re.search(r"([0-9]*D)", e.cigarstring),  "Make sure D's are merged to adjacent M"
                    reference_positions = e.get_reference_positions()  # Get genomic coordinates including 'D'

                try:
                    assigned_nucleotide = reference_positions[assignment]
                    # TODO: PERIODICTY GÖSTERMIYOR!
                except IndexError:  # if list index out of range
                    counter += 1
                    continue

                if e.reference_name not in assigned_positions:  # If the chromosome name is not present in the dictionary
                    assigned_positions[e.reference_name] = [assigned_nucleotide]
                else:  # Add the assigned position of the footprint into the relevant list in the dictionary
                    assigned_positions[e.reference_name].append(assigned_nucleotide)

        for chromosome in assigned_positions:  # When assignment ends for each entry in the chromosome dictionary
            assigned_positions[chromosome] = np.sort(
                np.array(assigned_positions[chromosome]))  # Make it np.array and sort

        if verbose and counter:
            print(f"{Col.WARNING}{counter} footprints were skipped because assignment index was out of range.{Col.ENDC}")

        return assigned_positions

    @staticmethod
    def footprint_counts_to_genes(footprint_genome_assignment_list, chromosome_gene_map, positions_gene_map, verbose=False):
        """
        Assigns the footprints to genes. Output shows how many footprints are assigned at each position of a gene.
        :param footprint_genome_assignment_list: List of outputs of footprint_assignment function
        :param chromosome_gene_map: Output of chromosomes_genes_matcher function
        :param positions_gene_map: Dictionary mapping gene_id to array of genomic coordinates, output of cds_ranges_to_array func.
        :param verbose: Boolean. Report or not report the progress.
        :return: Dictionary mapping gene_id to numpy array (with the same length of gene), which compose of integers.
        """
        max_possible = 2 ** 32  # Maximum assignment to a single nucleotide position. Since np.int32 is used.
        gene_footprint_assignment = dict()  # Initialize a dictionary to fill up

        # For each replicate do the same procedure
        for ind_assignment, footprint_genome_assignment in enumerate(footprint_genome_assignment_list):
            for ind, chr_name in enumerate(footprint_genome_assignment):  # For each chromosome
                if verbose:  # If 'verbose' activated, print the progress bar
                    progress_bar(ind, len(footprint_genome_assignment) - 1, suffix=f"    Chromosome: {chr_name}")
                # Get all genomic coordinates assigned to a footprint
                chr_pos_footprints = footprint_genome_assignment[chr_name]
                # Get the number of unique elements and the counts of these unique elements
                element_unique, element_count = np.unique(chr_pos_footprints, return_counts=True)
                # Create a dictionary from this information, mapping genomic position with the number of footprints assigned
                footprints_counts = dict(zip(element_unique, element_count))

                try:
                    gene_id_list = chromosome_gene_map[chr_name]  # For each gene in the same chromosome
                except KeyError:
                    continue
                    # Footprints, which are not mapped to a contig that has a gene on it, are excluded.

                for gene_id in gene_id_list:
                    pos_genes = positions_gene_map[gene_id]  # Get the genomic positions of the gene

                    # For each position in the gene, look at the footprints_counts, get the number if found, 0 otherwise.
                    temp_assignments = np.array([footprints_counts.get(i, 0) for i in pos_genes])

                    # Create an empty matrix to fill with raw counts, which originate from different replicates
                    if gene_id not in gene_footprint_assignment:
                        gene_footprint_assignment[gene_id] = np.zeros((len(footprint_genome_assignment_list), len(temp_assignments)),)
                    # Make sure the line is not filled previously
                    assert np.max(gene_footprint_assignment[gene_id][ind_assignment]) == 0, "Error in footprint_counts_to_genes: Multiple genes"
                    gene_footprint_assignment[gene_id][ind_assignment] = temp_assignments

        # Convert list of np.arrays to np.ndarray. Also set up the dtype as np.int32
        for gene_id in gene_footprint_assignment:
            assert np.max(gene_footprint_assignment[gene_id]) < max_possible, "Exceeded np.int32"
            gene_footprint_assignment[gene_id] = gene_footprint_assignment[gene_id].astype(np.int32)

        # There are some genes that are not annotated in any of the chromosomes or in the chromosomes in which there is no footprint.
        non_covered_genes = positions_gene_map.keys() - gene_footprint_assignment.keys()
        for gene_id in non_covered_genes:
            empty_matrix = np.zeros((len(footprint_genome_assignment_list), len(positions_gene_map[gene_id])), dtype=np.int32)
            assert gene_id not in gene_footprint_assignment, "Error in non_covered_genes"
            gene_footprint_assignment[gene_id] = empty_matrix

        return gene_footprint_assignment

    def calculate_rpm_genes(self, gene_id, average=True):
        replicates = self.total_assigned_gene[gene_id] / self.total_assigned * 1e6
        return np.mean(replicates) if average else replicates

    def calculate_rpkm_genes(self, gene_id, average=True):
        replicates = self.total_assigned_gene[gene_id] / self.total_assigned * 1e9 / self.gene_lengths[gene_id]
        return np.mean(replicates) if average else replicates

    def calculate_rpm_positions(self, gene_id, average=True):
        replicates = (self.gene_assignments[gene_id].T / self.total_assigned * 1e6).T
        return np.mean(replicates, axis=0) if average else replicates

    def calculate_confidence_interval_DEPRECATED(self, gene_id, value, confidence=0.95):
        assert value in ["lower", "upper", "mean"]
        rpm = self.calculate_rpm_positions(gene_id)
        mean = np.mean(rpm, axis=0)
        if value == "mean":
            return mean
        std_n_sqrt_conf = (np.std(rpm, axis=0) / np.sqrt(rpm.shape[0])) * confidence
        return mean - std_n_sqrt_conf if value == "lower" else mean + std_n_sqrt_conf

# todo: warnings.simplefilter('ignore', np.RankWarning)

class RiboSeqExperiment:

    def __init__(self, temp_repo_dir, sam_paths_translatome: list, sam_paths_experiment: list, name_experiment: str, assignment: int,
                 selection: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict, exclude_genes=[], verbose=True, recalculate=False):

        self.temp_repo_dir = temp_repo_dir
        self.sam_paths_translatome = sam_paths_translatome
        self.sam_paths_experiment = sam_paths_experiment
        self.name_experiment = name_experiment
        self.riboseq_assign_to = selection
        self.riboseq_assign_at = assignment
        self.exclude_genes = exclude_genes
        self.verbose = verbose
        self.recalculate = recalculate
        self.gene_list = sorted(gene_info_dictionary.keys())

        self.translatome = RiboSeqAssignment(self.sam_paths_translatome, self.temp_repo_dir, self.riboseq_assign_at,
                                             self.riboseq_assign_to, name_experiment + "_translatome",
                                             protein_genome_instance, gene_info_dictionary)
        self.experiment = RiboSeqAssignment(self.sam_paths_experiment, self.temp_repo_dir, self.riboseq_assign_at,
                                            self.riboseq_assign_to, name_experiment + "_experiment",
                                            protein_genome_instance, gene_info_dictionary)

        if exclude_genes:
            # Bu sayede kaydedilen assignment dosyaları tekrar tekrar kullanılabilir oluyor
            self.translatome.exclude_genes_calculate_stats(self.exclude_genes)
            self.experiment.exclude_genes_calculate_stats(self.exclude_genes)

    def normalized_rpm_positions(self, gene_id, smoothen=True):
        positions_experiment = self.experiment.calculate_rpm_positions(gene_id, average=True)
        positions_translatome = self.translatome.calculate_rpm_positions(gene_id, average=True)
        if smoothen:
            positions_experiment = smooth_array(positions_experiment, window_len=25, window="flat")
            positions_translatome = smooth_array(positions_translatome, window_len=25, window="flat")
        return positions_experiment / (positions_experiment + positions_translatome)


class RiboSeqSixtymers(RiboSeqExperiment):

    def __init__(self, temp_repo_dir, sam_paths_translatome: list, sam_paths_experiment: list, assignment: int,  # name_experiment: str,
                 selection: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict, exclude_genes=[], verbose=True, recalculate=False):
        super().__init__(temp_repo_dir, sam_paths_translatome, sam_paths_experiment, "sixtymers", assignment,
        selection, protein_genome_instance, gene_info_dictionary, exclude_genes=exclude_genes, verbose=verbose, recalculate=recalculate)

    def stalling_peaks_arpat(self, gene_id, mmc_threshold=0, normalized_peak_count_thr=0, get_top= 5):
        # Arbitrarily 5 for all from the Arpat paper.
        land = self.experiment.calculate_rpm_positions(gene_id, average=True)  # library size normalized disome peaks
        normalized_peak_count = np.sum(land)  # to test whether normalized peak count > 5
        mmc = np.mean(self.translatome.gene_assignments[gene_id])  # mean monosome count
        if mmc < mmc_threshold or normalized_peak_count < normalized_peak_count_thr:
            return np.nan
        else:
            normalized_rpm = land / mmc
            return normalized_rpm, normalized_rpm.argsort()[-get_top:][::-1]  # Arbitrarily 5

    def stalling_peaks_inecik_1(self):
        pass

    def see_examples(self, function, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            x, y = 3, 2
            fig, axes = plt.subplots(x, y, figsize=(y * 5, x), gridspec_kw={'hspace': 0, 'wspace': 0})
            random_gene_list = list()
            for i1 in range(axes.shape[0]):
                for i2 in range(axes.shape[1]):
                    gene_id = random.choice(self.gene_list)
                    arr, peaks = function(gene_id, *args, **kwargs)
                    total_exp = np.sum(self.experiment.total_assigned_gene[gene_id])
                    total_tra = np.sum(self.translatome.total_assigned_gene[gene_id])
                    random_gene_list.append(gene_id)
                    axes[i1][i2].plot(arr, alpha=0.75, color='gray')
                    axes[i1][i2].scatter(peaks, [arr[p] for p in peaks], color="blue", alpha=1, s=15)
                    axes[i1][i2].axes.get_xaxis().set_visible(False)
                    axes[i1][i2].axes.get_yaxis().set_visible(False)
                    axes[i1][i2].text(0, 0, f"{gene_id}: {total_exp} - {total_tra}", fontsize=6, transform=axes[i1][i2].transAxes)
            plt.tight_layout()
            plt.show()
        return random_gene_list


class RiboSeqSelective(RiboSeqExperiment):

    def __init__(self, temp_repo_dir, sam_paths_translatome: list, sam_paths_experiment: list, name_experiment: str, assignment: int,
                 selection: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict, exclude_genes=[], verbose=True, recalculate=False):
        super().__init__(temp_repo_dir, sam_paths_translatome, sam_paths_experiment, name_experiment, assignment,
        selection, protein_genome_instance, gene_info_dictionary, exclude_genes=exclude_genes, verbose=verbose, recalculate=recalculate)

    def calculate_binding_positions(self, normalized_rpm_threshold, min_gene_rpm_translatome, min_gene_rpkm_translatome):
        pass


class RiboSeqCoco(RiboSeqExperiment):

    def __init__(self, temp_repo_dir, sam_paths_monosome: list, sam_paths_disome: list, assignment: int,
                 selection: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict, exclude_genes=[], verbose=True, recalculate=False):
        super().__init__(temp_repo_dir, sam_paths_monosome, sam_paths_disome, "cocoassembly", assignment,
        selection, protein_genome_instance, gene_info_dictionary, exclude_genes=exclude_genes, verbose=verbose, recalculate=recalculate)

        self.output_file_name_fitting_calc = os.path.join(self.temp_repo_dir, f"riboseq_{self.name_experiment}_on_{self.riboseq_assign_to}_fitting_calculations.joblib")
        self.n_core = multiprocessing.cpu_count()

        try:
            assert os.access(self.output_file_name_fitting_calc, os.R_OK) and os.path.isfile(self.output_file_name_fitting_calc)
            if self.recalculate:
                print(
                    f"{Col.WARNING}Saved file is found at the path but 'recalculate=True': {self.output_file_name_fitting_calc}{Col.ENDC}.")
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
                    f"{Col.WARNING}There is at least one inconsistency between input parameters and loaded content. Recalculating...{Col.ENDC}")
                raise AssertionError
            self.best_model = loaded_content["best_model"]

        except (AssertionError, AttributeError, FileNotFoundError):
            print(f"{Col.HEADER}Sigmoid fitting are being calculated: {self.name_experiment}{Col.ENDC}")
            self.best_model = self.binomial_fitting(self.gene_list, n_core=self.n_core - 2)
            self.save_joblib()

    def save_joblib(self):
        # Write down the output dictionary and list as Joblib object for convenience in later uses.
        if self.verbose:
            print(f"{Col.HEADER}Calculations is being written to directory: {self.temp_repo_dir}{Col.ENDC}")
        to_dump = {"sam_paths_disome": self.sam_paths_experiment,
                   "sam_paths_monosome": self.sam_paths_translatome,
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
            print(f"{Col.HEADER}Fitting calculations found for {name_experiment} in path: {output_file_name}{Col.ENDC}")
        return joblib.load(output_file_name)

    def normalized_rpm_positions(self, *args, **kwargs):
        # Overrides parent method
        raise AssertionError("There is no translatome to calculate this.")

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
            executor.close()
            executor.join()
        result = [j for i in result for j in i]  # Flatten
        return dict(zip(gene_list, result))

    def binomial_fitting_single_core(self, gene_list):
        best_model = list()
        for ind, gene_id in enumerate(gene_list):
            if self.verbose:
                progress_bar(ind, len(gene_list) - 1)
            best_model.append(self.binomial_fitting_gene(gene_id))
        return best_model

    def create_bf(self, gene_id):  # BF = binomial fitting instance
        # below lines are identical to plot results
        disome_counts = np.sum(self.experiment.gene_assignments[gene_id], axis=0)  # mean problem çıkartıyor,
        monosome_counts = np.sum(self.translatome.gene_assignments[gene_id], axis=0)
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
            # todo: fonksiyon yerine açık açık yaz
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
            plt.show()
            if verbose:
                print(f"base::\nBic: {results['bic_scores']['base']}\nFun: {results['raw_fitting']['base'].fun}\n")
                print(f"ssig::\nBic: {results['bic_scores']['ssig']}\nFun: {results['raw_fitting']['ssig'].fun}\n")
                print(f"dsig::\nBic: {results['bic_scores']['dsig']}\nFun: {results['raw_fitting']['dsig'].fun}\n")
            print(f"Winner: {model_name}\nOnset: {onset}")

# TODO: benim coco datası 3-nucleotide periodicity gösteriyor mu
# TODO: KARAR: infrastructure'da neleri nasıl koyucaz burada


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


def biomart_mapping(temp_repo_dir, rscript, release=102):

    base_name = os.path.split(os.path.splitext(rscript)[0])[1]
    data_path = os.path.join(temp_repo_dir, f"{base_name}.txt")

    if not os.access(data_path, os.R_OK) or not os.path.isfile(data_path):
        print(f"{Col.HEADER}BiomaRt script is being run: {base_name}.{Col.ENDC}")
        r_installation = "RScript" if sys.platform == "darwin" else "Rscript"
        spr = subprocess.run(f"cd {temp_repo_dir}; {which(r_installation)} {rscript} {release} {data_path}", shell=True)
        assert spr.returncode == 0, f"Error: {rscript}"

    return pd.read_table(data_path)


def ensembl_release_object_creator(temp_repo_dir, release=102):

    def download_and_or_return(url_with_gz):
        basename_with_gz = os.path.split(url_with_gz)[1]
        basename_without_gz = os.path.splitext(os.path.split(url_with_gz)[1])[0]
        downloaded_file = os.path.join(temp_repo_dir, basename_with_gz)  # Downloaded gz path
        output_path = os.path.join(temp_repo_dir, basename_without_gz)  # Output path

        if not os.access(output_path, os.R_OK) or not os.path.isfile(output_path):  # Check if it already exist
            print(f"Downloading: {basename_with_gz}")
            spr1 = subprocess.run(f"cd {temp_repo_dir}; curl -L -R -O {url_with_gz}", shell=True)  # Download the file
            subprocess.run(f"cd {temp_repo_dir}; gzip -d -f {downloaded_file}", shell=True)  # Uncompress the file
            assert spr1.returncode == 0, f"Could not download: {url_with_gz}"

        return output_path

    gtf_url = f"ftp://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/" \
              f"Homo_sapiens.GRCh38.{release}.chr_patch_hapl_scaff.gtf.gz"
    transcript_fasta_url = f"ftp://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/" \
                           f"cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    protein_fasta_url = f"ftp://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/" \
                        f"pep/Homo_sapiens.GRCh38.pep.all.fa.gz"

    gtf_path = download_and_or_return(gtf_url)
    transcript_fasta_path = download_and_or_return(transcript_fasta_url)
    protein_fasta_path = download_and_or_return(protein_fasta_url)

    ero = pyensembl.Genome(
        reference_name="GRCh38",
        annotation_version=release,
        annotation_name=f"Homo_sapiens.GRCh38.{release}",
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
    strand = -1 if protein_genome_instance.db[transcript_list[0]][4] == "-" else 1  # all tx are the same
    # Ters olanları ters sırada veriyor!!!
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
    if not make_array:
        return cds_ranges
    else:
        strand = -1 if protein_genome_instance.db[best_transcript][4] == "-" else 1
        return np.concatenate([np.arange(i, j + strand, strand) for i, j in cds_ranges])


def ensembl_range_sum(ranges):
    """
    The function is to calculate the total CDS/exon length of a transcript.
    :param ranges: List of ranges. Obtained by coding_sequence_position_ranges or exon_intervals method of tx object
    :return: Integer. Total length of the CDS/exon
    """
    total = 0  # Initialize
    for i, j in ranges:  # For each CDS/exon of a transcript, get start and end positions
        # j always larger or equal to i
        total += (j - i + 1)  # Plus 1 is due to the fact that start and end is included.
    return total


def reduce_range_list(ranges):
    """
    Best way to return minimum number of ranges from a collection of ranges
    :param ranges: List of ranges (nested list), or list of tuples
    :return: List of lists
    """
    # Sorting the list based on the lower limit, and then based on size as tie breaker
    # Note: We need the bigger ranges to come before the smaller ranges
    rng = sorted(ranges, key=lambda pair: (pair[0], - pair[1]))
    reduced = []  # New set of ranges are stored here
    parent = rng[0]  # Use a parent range to decide if a range is inside other ranges
    # This will for sure be part of the solution, because it is the largest, leftmost range
    reduced.append(rng[0])  # Add it to the reduced list

    for x in rng:  # for each entry in ranges
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
        print(message)


class Col:
    HEADER = '\033[95m\033[1m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# End
