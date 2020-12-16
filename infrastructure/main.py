"""

"""


import os
import re
import subprocess
import sys
import warnings
import itertools

import joblib
import numpy as np
import pandas as pd
import pyensembl
import math
import pysam
from shutil import which
from Bio.Seq import translate
from Bio import pairwise2


infrastructure_script_path = __file__


class Infrastructre:

    def __init__(self, temp_repo_dir, exclude_genes=None, ensembl_release=102,
                 sam_translatome=None, sam_sixtymers=None,
                 include_gene3d=False, include_protein_genome=False):

        self.temp_repo_dir = temp_repo_dir
        self.exclude_genes = exclude_genes
        self.ensembl_release = ensembl_release

        self.gene_info_db_script = os.path.abspath(os.path.join(os.path.dirname(infrastructure_script_path), "gene_info_database.R"))
        gene_info_database = biomart_mapping(self.temp_repo_dir, self.gene_info_db_script, self.ensembl_release)
        gene_info_database = gene_info_database[gene_info_database["transcript_biotype"] == "protein_coding"].drop(columns=["transcript_biotype"])

        gene_list = np.unique(gene_info_database["ensembl_gene_id"].dropna())
        gene_list = np.setdiff1d(gene_list, exclude_genes)

        ero = ensembl_release_object_creator(self.temp_repo_dir, self.ensembl_release)

        self.gene_info = gene_class_dict_generate(self.temp_repo_dir, gene_list, gene_info_database, ero)

        if include_protein_genome:
            transcript_list = np.unique(gene_info_database["ensembl_transcript_id"].dropna())
            self.protein_genome = ProteinGenome(self.temp_repo_dir, transcript_list, ero)

        if include_gene3d:
            self.gene3d_db_script = os.path.abspath(os.path.join(os.path.dirname(infrastructure_script_path), "gene3d.R"))
            self.gene3d_database = EnsemblDomain(self.temp_repo_dir, self.gene3d_db_script, self.protein_genome, ero)

        if sam_translatome:
            self.translatome = RiboSeqAssignment(sam_translatome, self.temp_repo_dir, gene_list, -1, "translatome", ero)

        if sam_sixtymers:
            self.sixtymers = RiboSeqAssignment(sam_sixtymers, self.temp_repo_dir, gene_list, -1, "sixtymers", ero)

        # Uniprot
        # Conservation

    def calculate_stalling_peaks_arpat(self):
        try:
            pass
        except NameError:
            print(f"{Col.WARNING}Translatome and/or sixtymers are not defined.{Col.ENDC}")

    # annotation:
    # uniprot_annotations_ranges

    # conservation
    # conservation score -> ranges

    # get_riboseq, get_domain, get_conservation gibi function'lar yapılacakGE


def gene_class_dict_generate(temp_repo_dir, gene_list, gene_info_database, ensembl_release_object, overwrite=False, verbose=True):
    output_path = os.path.join(temp_repo_dir, "gene_info_database.joblib")
    if not os.access(output_path, os.R_OK) or not os.path.isfile(output_path) or overwrite:
        output = dict()
        if verbose:
            print(f"{Col.HEADER}Gene information dictionary are being created.{Col.ENDC}")
        for ind, gene_id in enumerate(gene_list):
            progress_bar(ind, len(gene_list) - 1, suffix=f"    {gene_id}")
            gi = gene_info_database[gene_info_database["ensembl_gene_id"] == gene_id]
            output[gene_id] = Gene(gi, ensembl_release_object)
        print(f"{Col.HEADER}Results are being written to directory: {temp_repo_dir}{Col.ENDC}")
        joblib.dump(output, output_path)
        if verbose:
            print(f"Done: {output_path}")
        return output
    else:
        if verbose:
            print(f"{Col.HEADER}Gene information dictionary is found in path: {output_path}{Col.ENDC}")
        return joblib.load(output_path)


class Gene:

    def __init__(self, one_gene_df, ensembl_release_object):

        temp_gene_id = one_gene_df["ensembl_gene_id"].unique()
        self.gene_id = temp_gene_id[0]
        self.gene_names = np.unique(np.concatenate((one_gene_df["external_gene_name"].dropna(), one_gene_df["external_synonym"].dropna())))
        temp_chromosome = np.unique(np.array(one_gene_df["chromosome_name"], dtype=str))
        self.chromosome = str(temp_chromosome[0])
        temp_strand = one_gene_df["strand"].unique()
        self.strand = "-" if temp_strand[0] == -1 else "+"
        temp_start = one_gene_df["start_position"].unique()
        self.start = int(temp_start[0])
        temp_end = one_gene_df["end_position"].unique()
        self.end = int(temp_end[0])
        assert all([len(i) == 1 for i in [temp_chromosome, temp_strand, temp_start, temp_end]])

        # Gene3d
        # pre_ensembl_peptide_ids = np.unique(one_gene_df["ensembl_peptide_id"])
        # pre_gene3d = gene3d_df[gene3d_df[""].isin(pre_ensembl_peptide_ids)].sort_values(
        #    by=["ensembl_peptide_id", "gene3d_start", "gene3d_end", "gene3d"], ignore_index=True)
        # Data
        # self.gene3d = pre_gene3d

        self.transcripts = self.prioritize_transcripts(one_gene_df)
        self.cds_ranges = gene_cds_ranges(ensembl_release_object, self.gene_id)

    @staticmethod
    def prioritize_transcripts(gene_df):
        transcripts = gene_df.drop(["ensembl_gene_id", "external_gene_name", "external_synonym",
                                    "chromosome_name", "start_position", "end_position", "strand"], axis=1).drop_duplicates()
        transcripts["transcript_appris"] = transcripts["transcript_appris"].replace(
            ["alternative1", "alternative2"], ["renamed_alternative1", "renamed_alternative2"])
        # apris'de üçbeş
        assert len(np.unique(transcripts["transcript_mane_select"].dropna())) <= 1
        # mane_clinical var mı yok my
        # basic or not
        # tsl sıralı zaten
        transcripts.sort_values(by=["transcript_mane_select", "transcript_appris", "transcript_gencode_basic",
                                    "transcript_tsl", "external_transcript_name"], inplace=True, ignore_index=True)
        return transcripts


class ProteinGenome:

    def __init__(self, temp_repo_dir, transcript_list, ensembl_release_object, verbose=True, recalculate=False):
        self.transcript_list = np.sort(transcript_list)
        self.ensembl_release_object = ensembl_release_object
        self.temp_repo_dir = temp_repo_dir
        self.verbose = verbose
        self.file_name = "protein_genome_instance.joblib"
        self.output_file_name = os.path.join(self.temp_repo_dir, self.file_name)
        self.recalculate = recalculate

        try:
            assert os.access(self.output_file_name, os.R_OK) and os.path.isfile(self.output_file_name)
            if self.recalculate:
                print(f"{Col.WARNING}Saved file is found at the path but 'recalculate=True': {self.output_file_name}{Col.ENDC}.")
                raise AssertionError
            loaded_content = self.load_joblib(self.output_file_name, self.verbose)
            consistent_with = all([
                (self.transcript_list == loaded_content.transcript_list).all(),
                self.ensembl_release_object == loaded_content.ensembl_release_object,
                self.temp_repo_dir == loaded_content.temp_repo_dir,
                self.file_name == loaded_content.file_name,
                self.output_file_name == loaded_content.output_file_name,
            ])
            if not consistent_with:
                print(f"{Col.WARNING}There is at least one inconsistency between input parameters and loaded content. Recalculating...{Col.ENDC}")
                raise AssertionError
            self.transcripts = loaded_content.transcripts
        except (AssertionError, FileNotFoundError):
            print(f"{Col.HEADER}Protein genome mapping are being calculated.{Col.ENDC}")
            self.calculate_transcript_mapping()
            self.save_joblib()

    @staticmethod
    def consistent_coding_ranges(transcript_object):

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

        coding_position_ranges = [list(i) for i in transcript_object.coding_sequence_position_ranges]
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
            query_transcript = query_transcript + "N" * modify3 if modify3 >=0 else query_transcript[:modify3]

            assert ensembl_range_sum(coding_position_ranges) == len(transcript_object.protein_sequence) * 3
            assert len(query_transcript) == len(transcript_object.protein_sequence) * 3

            return [[i, j] if transcript_object.strand == "+" else [j, i] for i, j in coding_position_ranges], \
                   transcript_object.protein_sequence, query_transcript, transcript_object.contig, transcript_object.strand

    def calculate_transcript_mapping(self):
        assert len(np.unique(self.transcript_list)) == len(self.transcript_list)
        output = dict()
        for ind, transcript_id in enumerate(self.transcript_list):
            if self.verbose:
                progress_bar(ind, len(self.transcript_list) - 1)
            transcript_object = self.ensembl_release_object.transcript_by_id(transcript_id)
            output[transcript_id] = ProteinGenome.consistent_coding_ranges(transcript_object)
        self.transcripts = output

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
        transcript = self.transcripts[transcript_id]

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
        for ind , (cdr_start, cdr_end) in enumerate(cdr_relative):
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
                self.temp_repo_dir == loaded_content.temp_repo_dir,
                self.rscript == loaded_content.rscript,
                self.base_name == loaded_content.base_name,
                self.data_path == loaded_content.data_path,
                self.ensembl_release == loaded_content.ensembl_release,
                self.output_file_name == loaded_content.output_file_name,
            ])
            if not consistent_with:
                print(f"{Col.WARNING}There is at least one inconsistency between input parameters and loaded content. Recalculating...{Col.ENDC}")
                raise AssertionError
            self.df = loaded_content.df
            self.columns = loaded_content.columns
        except (AssertionError, FileNotFoundError):
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
        print(self.columns)
        coordinates_ranges = list()
        coordinates_contig = list()
        for ind, (protein_id, domain_start, domain_end) in enumerate(zip(protein_ids, domain_starts, domain_ends)):
            if self.verbose:
                progress_bar(ind, len(protein_ids) - 1)
            transcript_id = ensembl_release_object.transcript_id_of_protein_id(protein_id)
            coordinates_ranges.append(protein_genome_instance.protein2genome(transcript_id, domain_start, domain_end))
            coordinates_contig.append(protein_genome_instance.transcripts[transcript_id][3])
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


# TODO: LIST
# uzaktan erişimin mümkün olup olmadığını sor: IT guy'a
# bol bol açıklama ekle. bütün koda (protein_sequence asıl: nedenini açıkla)
# uniprot ve conservation'u tatilde yap
# gene class'ını bitir.
# co-co site'ı hesaplayan fitting şeyi yaz: coco_site'ı eklemlendir.
#
# TODO: TOMORROW
# çocuğun analizine başla


class RiboSeqAssignment:

    def __init__(self, sam_paths, temp_repo_dir, gene_list, assignment, riboseq_group, ensembl_release_object, verbose=True, recalculate=False):
        # :param riboseq_group: String to identify the RiboSeq experiment annotation; like "translatome" or "60mers"
        self.sam_paths = np.sort(sam_paths)
        self.temp_repo_dir = temp_repo_dir
        self.output_file_name = os.path.join(temp_repo_dir, f"riboseq_on_genes_{riboseq_group}.joblib")
        self.gene_list = np.sort(gene_list)
        self.assignment = assignment
        self.verbose = verbose
        self.riboseq_group = riboseq_group  # Define riboseq_group variable for a function.
        self.ensembl_release_object = ensembl_release_object
        self.recalculate = recalculate

        try:
            assert os.access(self.output_file_name, os.R_OK) and os.path.isfile(self.output_file_name)
            if self.recalculate:
                print(f"{Col.WARNING}Saved file is found at the path but 'recalculate=True': {self.output_file_name}{Col.ENDC}.")
                raise AssertionError
            loaded_content = self.load_joblib(self.riboseq_group, self.output_file_name, self.verbose)
            consistent_with = all([
                (self.sam_paths == loaded_content.sam_paths).all(),
                self.temp_repo_dir == loaded_content.temp_repo_dir,
                self.output_file_name == loaded_content.output_file_name,
                self.riboseq_group == loaded_content.riboseq_group,
                (self.gene_list == loaded_content.gene_list).all(),
                self.assignment == loaded_content.assignment,
                self.ensembl_release_object == loaded_content.ensembl_release_object
            ])
            if not consistent_with:
                print(f"{Col.WARNING}There is at least one inconsistency between input parameters and loaded content. Recalculating...{Col.ENDC}")
                raise AssertionError
            self.gene_assignments = loaded_content.gene_assignments
            self.gene_lengths = loaded_content.gene_lengths
            self.total_assigned_gene = loaded_content.total_assigned_gene
            self.total_assigned = loaded_content.total_assigned

        except (AssertionError, FileNotFoundError):
            print(f"{Col.HEADER}Gene assignments are being calculated: {self.riboseq_group}{Col.ENDC}")
            self.calculate_gene_assignments()
            self.gene_lengths = np.array([self.gene_assignments[i].shape[1] if i in self.gene_assignments else np.nan for i in self.gene_list])
            self.total_assigned_gene = np.array([np.sum(self.gene_assignments[i], axis=1) if i in self.gene_assignments else np.zeros(len(self.sam_paths)) for i in self.gene_list])
            self.total_assigned = np.sum(self.total_assigned_gene)
            self.save_joblib()

    def get_assigned_read_mapped_read_all(self):
        mapped_read = 0
        for sam_path in self.sam_paths:
            with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open SAM file with pysam library
                mapped_read += sam_handle.count()  # Call count method to find out total read count.

        return self.total_assigned, mapped_read

    def get_assigned_read_mapped_read_individual(self):
        mapped_read = list()
        for sam_index in self.sam_paths:
            with pysam.AlignmentFile(sam_index, "r") as sam_handle:  # Open SAM file with pysam library
                mapped_read.append(sam_handle.count())  # Call count method to find out total read count.
        assigned_read = list(np.sum(self.total_assigned_gene, axis=0))
        return {l: [a, m] for a, m, l in zip(assigned_read, mapped_read, list(self.sam_paths))}

    def calculate_gene_assignments(self):
        if self.verbose:
            print(f"{Col.HEADER}Genes are allocated to chromosomes.{Col.ENDC}")
        # Get the mapping from chromosomes to list of genes which they carry
        chromosome_gene = self.chromosomes_genes_matcher(self.ensembl_release_object, self.gene_list)
        if self.verbose:
            print(f"{Col.HEADER}CDS ranges of genes are being calculated.{Col.ENDC}")
        # Get array of genomic positions for all genes, create a dictionary out of it.
        positions_gene = {i: gene_cds_array(self.ensembl_release_object, i) for i in self.gene_list}

        if self.verbose:
            print(f"{Col.HEADER}Footprint are being assigned to genomic coordinates.{Col.ENDC}")
        # Assign footprints to genomic positions
        footprint_genome_assignment_list = [self.footprint_assignment(sam_path, assignment=self.assignment, verbose=self.verbose) for sam_path in self.sam_paths]
        # Remove chromosomes which does not contain footprint assignment, or which does not contain genes
        for i in range(len(footprint_genome_assignment_list)):
            chromosome_gene, footprint_genome_assignment_list[i] = sync_dictionaries(chromosome_gene, footprint_genome_assignment_list[i])

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
    def chromosomes_genes_matcher(ensembl_release_object, gene_list):
        """
        This script is to assign genes into chromosomes.
        :param ensembl_release_object: Created by ensembl_release_object_creator() function
        :param gene_list: List of genes to assign chromosomes
        :return: Dictionary of assigned genes.
        """
        chromosomes = dict()  # Initialize a dictionary to fill up
        for gene in gene_list:  # For each gene in the gene_list
            gene_object = ensembl_release_object.gene_by_id(gene)  # Create an instance of Ensembl Gene object
            # Find the chromosome name (contig), add the gene into the relevant list
            chromosomes[gene_object.contig] = chromosomes.get(gene_object.contig, []) + [gene_object.id]
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
                    continue  # Genes which are not among gene_list are skipped here as well.
                    # For example I saw ENSG00000278457, gives disome footprint but it is not among gene_list
                    # chromosome_gene_map is based on gene_list only.

                for gene_id in gene_id_list:
                    pos_genes = positions_gene_map[gene_id]  # Get the genomic positions of the gene
                    # For each position in the gene, look at the footprints_counts, get the number if found, 0 otherwise.
                    gene_footprint_assignment[gene_id] = gene_footprint_assignment.get(gene_id, []) + [np.array([footprints_counts.get(i, 0) for i in pos_genes])]

        # Convert list of np.arrays to np.ndarray. Also set up the dtype as np.int32
        for gene_id in gene_footprint_assignment:
            assert all([np.max(i) < max_possible for i in gene_footprint_assignment[gene_id]]), "Exceeded np.int32"
            gene_footprint_assignment[gene_id] = np.stack(gene_footprint_assignment[gene_id], axis=0).astype(np.int32)

        return gene_footprint_assignment

# todo: dictionary yerine 3d array yap
# todo: assigned mapped ayrımı yap isimlendirmede

    def calculate_rpm_genes(self):
        return self.total_assigned_gene / self.total_assigned * 10**6

    def calculate_rpkm_genes(self):
        return self.total_assigned_gene / self.total_assigned * 10**9 / np.stack([self.gene_lengths] * len(self.sam_paths), axis=1)

    def calculate_rpm_positions(self, gene_id):
        return self.gene_assignments[gene_id] / self.total_assigned * 10**6

    def calculate_confidence_interval(self, gene_id, value, confidence=0.95):
        assert value in ["lower", "upper", "mean"]
        rpm = self.calculate_rpm_positions(gene_id)
        mean = np.mean(rpm, axis=0)
        if value == "mean":
            return mean
        std_n_sqrt_conf = (np.std(rpm, axis=0) / np.sqrt(rpm.shape[0])) * confidence
        return mean - std_n_sqrt_conf if value == "lower" else mean + std_n_sqrt_conf


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
        spr = subprocess.run(f"cd {temp_repo_dir}; {which('RScript')} {rscript} {release} {data_path}", shell=True)
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


def gene_cds_ranges(ensembl_release_object, gene_id):
    cds_ranges = ensembl_release_object.db.query(select_column_names=["start", "end"], feature="CDS", distinct=True,
                                                 filter_column="gene_id", filter_value=gene_id)
    return reduce_range_list(cds_ranges)  # Warning: it also sorts!!!


def gene_cds_array(ensembl_release_object, gene_id):
    cds_ranges = gene_cds_ranges(ensembl_release_object, gene_id)
    # +1 → both ends should be included.
    return np.concatenate([np.arange(i, j + 1) for i, j in cds_ranges])


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


def progress_bar(iteration, total, prefix='Progress:', suffix='', decimals=1, bar_length=20):
    """
    This function should be called inside of loop, gives the loop's progress.
    :param iteration: It is integer. It is current iteration.
    :param total: It is integer. It is total iteration.
    :param prefix: It is string. It will be placed before progress bar.
    :param suffix: It is string. It will be placed after progress bar.
    :param decimals: It is integer. It is number of decimals in percent complete.
    :param bar_length: It is integer. It is character length of bar.
    :return: It is void function. Nothing is returned.
    """
    filled_length = int(round(bar_length * iteration / float(total)))
    percents = round(100.00 * (iteration / float(total)), decimals)
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()


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


if __name__ == '__main__':
    sam_paths_translatome = [
        "/Users/kemalinecik/Documents/dropbox/TT1.sam",
        "/Users/kemalinecik/Documents/dropbox/TT2.sam",
    ]
    sam_paths_sixtymers = [
        "/Users/kemalinecik/Documents/dropbox/Rep1.sam",
        "/Users/kemalinecik/Documents/dropbox/Rep2.sam",
        "/Users/kemalinecik/Documents/dropbox/NoPK.sam"
    ]
    exclude_genes = np.array(["ENSG00000160789"])
    temp_repo_dir = "/Users/kemalinecik/Documents/dropbox"
    Infrastructre(temp_repo_dir, exclude_genes)


# End
