"""

"""


import sys
import os
import re
import subprocess
import sys
import warnings
import itertools

import joblib
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing
from scipy.optimize import basinhopping
from statsmodels.stats.proportion import proportion_confint
import pyensembl
import math
import pysam
from shutil import which
from Bio.Seq import translate
from Bio import pairwise2
import matplotlib.pyplot as plt


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

        transcript_list = np.unique(gene_info_database["ensembl_transcript_id"].dropna())
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


def gene_class_dict_generate(temp_repo_dir, gene_list, gene_info_database, gene_info_names, gene_info_uniprot,
                             overwrite=False, verbose=True):

    output_path = os.path.join(temp_repo_dir, "gene_info_database.joblib")
    if not os.access(output_path, os.R_OK) or not os.path.isfile(output_path) or overwrite:
        output = dict()
        if verbose:
            print(f"{Col.HEADER}Gene information dictionary are being created.{Col.ENDC}")
        for ind, gene_id in enumerate(gene_list):
            progress_bar(ind, len(gene_list) - 1, suffix=f"    {gene_id}")
            gi = gene_info_database[gene_info_database["ensembl_gene_id"] == gene_id]
            gi_names = gene_info_names[gene_info_names["ensembl_gene_id"] == gene_id]
            gi_uniprot = gene_info_uniprot[gene_info_uniprot["ensembl_gene_id"] == gene_id]
            output[gene_id] = Gene(gi, gi_names, gi_uniprot)
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

    def __init__(self, one_gene_df, one_gene_names, one_gene_uniprot):

        temp_gene_id = one_gene_df["ensembl_gene_id"].unique()
        self.gene_id = temp_gene_id[0]
        self.gene_names = np.unique(np.concatenate((one_gene_names["external_gene_name"].dropna(), one_gene_names["external_synonym"].dropna())))
        self.uniprot_ids = np.unique(np.concatenate((one_gene_uniprot["uniprotswissprot"].dropna(), one_gene_uniprot["uniprotsptrembl"].dropna())))
        temp_chromosome = np.unique(np.array(one_gene_df["chromosome_name"], dtype=str))
        self.chromosome = str(temp_chromosome[0])
        temp_strand = one_gene_df["strand"].unique()
        self.strand = "-" if temp_strand[0] == -1 else "+"
        temp_start = one_gene_df["start_position"].unique()
        self.start = int(temp_start[0])
        temp_end = one_gene_df["end_position"].unique()
        self.end = int(temp_end[0])
        assert all([len(i) == 1 for i in [temp_chromosome, temp_strand, temp_start, temp_end]])
        self.transcripts = self.prioritize_transcripts(one_gene_df)

    def prioritize_transcripts(self, gene_df):
        transcripts = gene_df.drop(["ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"], axis=1).drop_duplicates()
        transcripts["transcript_appris"] = transcripts["transcript_appris"].replace(
            ["alternative1", "alternative2"], ["renamed_alternative1", "renamed_alternative2"])
        # apris'de üçbeş
        assert len(np.unique(transcripts["transcript_mane_select"].dropna())) <= 1
        # mane_clinical var mı yok my
        # basic or not
        # tsl sıralı zaten

        transcripts.sort_values(by=["transcript_mane_select", "transcript_appris", "transcript_gencode_basic",
                                    "transcript_tsl", "external_transcript_name"], inplace=True, ignore_index=True)
        # apris uzun olanı seçiyor zaten
        return transcripts


class ProteinGenome:

    def __init__(self, temp_repo_dir, transcript_list, ensembl_release_object, verbose=True, recalculate=False):
        self.transcript_list = sorted(transcript_list)
        self.ensembl_release = ensembl_release_object.annotation_version
        self.temp_repo_dir = temp_repo_dir
        self.verbose = verbose
        self.output_file_name = os.path.join(self.temp_repo_dir, "protein_genome_instance.joblib")
        self.recalculate = recalculate

        try:
            assert os.access(self.output_file_name, os.R_OK) and os.path.isfile(self.output_file_name)
            if self.recalculate:
                print(f"{Col.WARNING}Saved file is found at the path but 'recalculate=True': {self.output_file_name}{Col.ENDC}.")
                raise AssertionError
            loaded_content = self.load_joblib(self.output_file_name, self.verbose)
            consistent_with = all([
                self.transcript_list == loaded_content.transcript_list,
                self.ensembl_release == loaded_content.ensembl_release,
                os.path.basename(self.output_file_name) == os.path.basename(loaded_content.output_file_name),
            ])
            if not consistent_with:
                print(f"{Col.WARNING}There is at least one inconsistency between input parameters and loaded content. Recalculating...{Col.ENDC}")
                raise AssertionError
            self.db = loaded_content.db
        except (AssertionError, AttributeError, FileNotFoundError):
            print(f"{Col.HEADER}Protein genome mapping are being calculated.{Col.ENDC}")
            self.calculate_transcript_mapping(ensembl_release_object)
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

        try:
            coding_position_ranges = [list(i) for i in transcript_object.coding_sequence_position_ranges]
        except ValueError: # ValueError: Transcript does not contain feature CDS
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
            query_transcript = query_transcript + "N" * modify3 if modify3 >=0 else query_transcript[:modify3]

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
        self.db = output

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
                 protein_genome_instance:ProteinGenome, gene_info_dictionary:dict, verbose=True, recalculate=False):
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
            self.gene_lengths = np.array([int(self.gene_assignments[i].shape[1]) for i in self.gene_list], dtype=int)
            self.exclude_genes_calculate_stats(self.exclude_gene_list)
            self.save_joblib()

    def exclude_genes_calculate_stats(self, exclude_gene_list:list):
        for gene_id in exclude_gene_list:
            target_shape = self.gene_assignments[gene_id].shape
            empty_matrix = np.empty(target_shape)
            empty_matrix.fill(np.nan)
            self.gene_assignments[gene_id] = empty_matrix

        self.exclude_gene_list.extend(exclude_gene_list)
        self.total_assigned_gene = np.array([np.sum(self.gene_assignments[i], axis=1) for i in self.gene_list], dtype=int)
        self.total_assigned = int(np.sum(self.total_assigned_gene))

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
                        gene_footprint_assignment[gene_id] = np.zeros((len(footprint_genome_assignment_list),len(temp_assignments)),)
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

    def calculate_rpm_genes(self):
        return self.total_assigned_gene / self.total_assigned * 10**6

    def calculate_rpkm_genes(self):
        return self.total_assigned_gene / self.total_assigned * 10**9 / np.stack([self.gene_lengths] * len(self.sam_paths), axis=1)

    def calculate_rpm_positions(self, gene_id):
        return self.gene_assignments.get(gene_id, np.zeros(self.gene_lengths[np.where(self.gene_list==gene_id)])) / self.total_assigned * 10**6

    def calculate_confidence_interval_DEPRECATED(self, gene_id, value, confidence=0.95):
        assert value in ["lower", "upper", "mean"]
        rpm = self.calculate_rpm_positions(gene_id)
        mean = np.mean(rpm, axis=0)
        if value == "mean":
            return mean
        std_n_sqrt_conf = (np.std(rpm, axis=0) / np.sqrt(rpm.shape[0])) * confidence
        return mean - std_n_sqrt_conf if value == "lower" else mean + std_n_sqrt_conf


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

    def normalized_rpm_for_positions(self, min_gene_rpm_translatome=0, min_gene_rpkm_translatome=0):
        output = dict()
        translatome_rpm_mean = np.mean(self.translatome.calculate_rpm_genes(), axis=1)
        translatome_rpkm_mean = np.mean(self.translatome.calculate_rpkm_genes(), axis=1)
        for ind, gene_id in enumerate(self.gene_list):  # gene_info_dictionary.keys()?
            position_rpm_mean = np.mean(self.experiment.calculate_rpm_positions(gene_id), axis=0)
            gene_rpm_mean = translatome_rpm_mean[ind]
            normalized_rpm = position_rpm_mean / gene_rpm_mean
            if gene_rpm_mean < min_gene_rpm_translatome or translatome_rpkm_mean[ind] < min_gene_rpkm_translatome:
                output[gene_id] = np.nan
            else:
                output[gene_id] = normalized_rpm
        return output

class RiboSeqSixtymers(RiboSeqExperiment):

    def __init__(self, temp_repo_dir, sam_paths_translatome: list, sam_paths_experiment: list, assignment: int,  # name_experiment: str,
                 selection: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict, exclude_genes=[], verbose=True, recalculate=False):
        super().__init__(temp_repo_dir, sam_paths_translatome, sam_paths_experiment, "sixtymers", assignment,
        selection, protein_genome_instance, gene_info_dictionary, exclude_genes = exclude_genes, verbose = verbose, recalculate = recalculate)


    def calculate_stalling_peaks_arpat(self):
        output = dict()
        translatome_rpm_mean = np.mean(self.translatome.calculate_rpm_genes(), axis=1)
        for ind, gene_id in enumerate(self.translatome.gene_list):  # gene_info_dictionary.keys()?
            position_rpm_mean = np.mean(self.experiment.calculate_rpm_positions(gene_id), axis=0)
            normalized_peak_count = np.sum(position_rpm_mean != 0)
            gene_rpm_mean = translatome_rpm_mean[ind]
            if gene_rpm_mean < 5 or normalized_peak_count < 5:  # Arbitrarily 5
                output[gene_id] = np.nan
            else:
                normalized_rpm = position_rpm_mean / gene_rpm_mean
                output[gene_id] = normalized_rpm.argsort()[-5:][::-1] # Arbitrarily 5
                # todo: genome position'larını return etsin
        return output

    def calculate_stalling_peaks_inecik(self):
        pass  # todo


class RiboSeqSelective(RiboSeqExperiment):

    def __init__(self, temp_repo_dir, sam_paths_translatome: list, sam_paths_experiment: list, name_experiment: str, assignment: int,
                 selection: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict, exclude_genes=[], verbose=True, recalculate=False):
        super().__init__(temp_repo_dir, sam_paths_translatome, sam_paths_experiment, name_experiment, assignment,
        selection, protein_genome_instance, gene_info_dictionary, exclude_genes = exclude_genes, verbose = verbose, recalculate = recalculate)

    def calculate_binding_positions(self, normalized_rpm_threshold, min_gene_rpm_translatome, min_gene_rpkm_translatome):
        normalized_rpm = self.normalized_rpm_for_positions(min_gene_rpm_translatome, min_gene_rpkm_translatome)
        output = dict()
        for ind, gene_id in enumerate(normalized_rpm.keys()):
            gene_normalized_rpm = output[gene_id]
            if np.isnan(gene_normalized_rpm):
                output[gene_id] = np.nan
            else:
                output[gene_id] = np.where(normalized_rpm[gene_id] > normalized_rpm_threshold)
            # todo: genome position'larını return etsin
        return output


class RiboSeqCoco(RiboSeqExperiment):

    def __init__(self, temp_repo_dir, sam_paths_monosome: list, sam_paths_disome: list, assignment: int,
                 selection: str, protein_genome_instance: ProteinGenome, gene_info_dictionary: dict, exclude_genes=[], verbose=True, recalculate=False):
        super().__init__(temp_repo_dir, sam_paths_monosome, sam_paths_disome, "cocoassembly", assignment,
        selection, protein_genome_instance, gene_info_dictionary, exclude_genes = exclude_genes, verbose = verbose, recalculate = recalculate)

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
            self.binomial_fitting(self.gene_list, n_core=self.n_core - 2) # todo:CORRECT LINE
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

    def normalized_rpm_for_positions(self, *args, **kwargs):
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
        self.best_model = dict(zip(gene_list, result))
        # Sıralı geliyor mu kontrol et bi

    def binomial_fitting_single_core(self, gene_list):
        best_model = list()
        for ind, gene_id in enumerate(gene_list):
            if self.verbose:
                progress_bar(ind, len(gene_list) - 1)
            best_model.append(self.binomial_fitting_gene(gene_id))
        return best_model

    def binomial_fitting_gene(self, gene_id):
        disome_counts = np.sum(self.experiment.gene_assignments[gene_id], axis=0)  # mean problem çıkartıyor,
        monosome_counts = np.sum(self.translatome.gene_assignments[gene_id], axis=0)
        x_data = np.arange(1, len(disome_counts) + 1)
        fitter_instance = BinomialFitting(x_data, disome_counts, monosome_counts)
        model_type, parameters = fitter_instance.select_correct_model()
        return model_type, parameters, len(disome_counts)

    def calculate_curve(self, gene_id):
        model_type, parameters, x_length = self.best_model[gene_id]
        x_data = np.arange(1, x_length + 1)
        if model_type == "base":
            y_predicted = BinomialFitting.model_base(x_data, *parameters)
        elif model_type == "ssig":
            y_predicted = BinomialFitting.model_ssig(x_data, *parameters)
        elif model_type == "dsig":
            y_predicted = BinomialFitting.model_dsig(x_data, *parameters)
        else:
            raise ValueError("Unknown model in calculate_curve function.")
        return y_predicted

    def calculate_onset(self,gene_id):
        model_type, parameters, x_data = self.best_model[gene_id]
        if model_type == "base":
            return np.nan
        elif model_type == "ssig":
            return parameters[3]  # i_mid
        elif model_type == "dsig":
            return parameters[5]  # i_mid
        else:
            raise ValueError("Unknown model in calculate_onset function.")

    def plot_result(self, gene_id):
        disome_counts = np.sum(self.experiment.gene_assignments[gene_id], axis=0)  # mean problem çıkartıyor,
        monosome_counts = np.sum(self.translatome.gene_assignments[gene_id], axis=0)
        x_data = np.arange(1, len(disome_counts) + 1)
        fitter = BinomialFitting(x_data, disome_counts, monosome_counts)
        fitter.plot_result()
        return fitter


# TODO: benim fitting aa değil nükleotid alıyor bu sorun olur mu, makalede ne diyor buna dair?
# TODO: benim coco datası 3-nucleotide periodicity gösteriyor mu
# TODO: onset'ler korrele mi makaledekiyle?

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
        return p

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

    def neg_log_likelihood_base(self,param):
        # The same calculation at model_base(), but for all x_data
        y_predicted = self.model_base(self.x_data, *param)
        negative_log_likelihood = -np.sum(stats.binom.logpmf(k=self.disome, n=self.n_trial, p=y_predicted))
        return negative_log_likelihood

    def neg_log_likelihood_ssig(self,param):
        # The same calculation at model_ssig(), but for all x_data
        y_predicted = self.model_ssig(self.x_data, *param)
        negative_log_likelihood = -np.sum(stats.binom.logpmf(k=self.disome, n=self.n_trial, p=y_predicted))
        return negative_log_likelihood

    def neg_log_likelihood_dsig(self, param):
        # The same calculation at model_dsig(), but for all x_data
        y_predicted = self.model_dsig(self.x_data, *param)
        negative_log_likelihood = -np.sum(stats.binom.logpmf(k=self.disome, n=self.n_trial, p=y_predicted))
        return negative_log_likelihood

    def minimize_base(self, niter=100, seed=1):
        x0 = np.array([0.5])
        bounds = ((0,1),)
        minimizer_kwargs = {"method": "SLSQP", "bounds": bounds}
        return basinhopping(func=self.neg_log_likelihood_base, x0=x0, minimizer_kwargs=minimizer_kwargs,
                            seed=seed, niter=niter)

    def minimize_ssig(self, niter=100, seed=1):
        x0 = np.array([0.25, 0.65, 0.25, int(len(self.x_data)/2)])
        bounds = ((0, 1), (0, 1), (0, 0.5), (1, len(self.x_data)))
        minimizer_kwargs = {"method": "SLSQP", "bounds": bounds}
        return basinhopping(func=self.neg_log_likelihood_ssig, x0=x0, minimizer_kwargs=minimizer_kwargs,
                            seed=seed, niter=niter)

    def minimize_dsig(self, niter=100, seed=1):
        x0 = np.array([0.25, 0.75, 0.5, 0.15, -0.3, int(len(self.x_data) / 2), int(len(self.x_data) / 2)])
        bounds = ((0, 1), (0, 1), (0, 1), (0, 0.5), (-0.5, 0), (1, len(self.x_data)), (1, len(self.x_data)))
        minimizer_kwargs = {"method": "SLSQP", "bounds": bounds}
        return basinhopping(func=self.neg_log_likelihood_dsig, x0=x0, minimizer_kwargs=minimizer_kwargs,
                            seed=seed, niter=niter)

    def minimize_models(self):
        return self.minimize_base(), self.minimize_ssig(), self.minimize_dsig()

    def calculate_BIC(self, minimized_result):  # Bayesian Information Criterion
        # BIC = -2 * LL + log(N) * k
        # Where log() has the base-e called the natural logarithm, LL is the log-likelihood of the model, N is the number of examples in the training dataset, and k is the number of parameters in the model.
        # The score as defined above is minimized, e.g. the model with the lowest BIC is selected.
        k = len(minimized_result.x)
        LL = -minimized_result.fun
        N = len(self.x_data)
        return -2 * LL + np.log(N) * k

    def select_correct_model(self):
        minimized_results = self.minimize_models()
        bic = [self.calculate_BIC(i) for i in minimized_results]
        winner = sorted(zip(bic, ["base", "ssig", "dsig"], minimized_results))[0]
        return winner[1], winner[2].x

    def plot_result(self):
        base, ssig, dsig = self.minimize_models()
        plt.plot(self.model_base(self.x_data, *base.x) + np.zeros(len(self.x_data)), color='black')
        plt.plot(self.model_ssig(self.x_data, *ssig.x), color='blue')
        plt.plot(self.model_dsig(self.x_data, *dsig.x), color='red')
        plt.scatter(self.x_data, self.disome / (self.disome + self.monosome),alpha=0.25, color='gray')
        plt.show()
        print(f"Base ->\nBic:{self.calculate_BIC(base)}\nFun: {base.fun}\n")
        print(f"Ssig ->\nBic:{self.calculate_BIC(ssig)}\nFun: {ssig.fun}\n")
        print(f"Dsig ->\nBic:{self.calculate_BIC(dsig)}\nFun: {dsig.fun}\n")

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
    if strand == -1:
        cds_all_ranges = [[j, i] for i, j in reversed(cds_all_ranges)]  # reverse if at '-'
    if not make_array:
        return cds_all_ranges
    else:
        return np.concatenate([np.arange(i, j + strand, strand) for i, j in cds_all_ranges])


def best_transcript_cds(protein_genome_instance, gene_info_dictionary, gene_id, make_array=False):
    best_transcript = gene_info_dictionary[gene_id].transcripts.iloc[0][0]  # At least 1 transcript exists
    cds_ranges = protein_genome_instance.db[best_transcript][0]
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


def gtf_cds_parser(gtf_path, verbose=False):
    """
    This is to convert gtf file into list of dictionary elements.
    :param verbose: Boolean. Report or not report the progress.
    :param gtf_path: String. Path of the gtf or gtf.gz file
    :return: List of dictionaries for each entry
    """
    # Column names for GTF file. Source: https://www.ensembl.org/info/website/upload/gff.html
    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame"]
    assert os.path.splitext(gtf_path)[1] == ".gtf"
    with open(gtf_path, "r") as gtf_handle:  # Open and read the file
        gtf = gtf_handle.readlines()
    # Remove the titles, strip the trailing white spaces, split the line into cells
    gtf = [i.strip().split('\t') for i in gtf if not i.startswith('#')]
    output = list()
    filter_columns = ["protein_id", "transcript_id", "gene_name", "start", "end", "strand", "frame", "seqname", "exon_number"]
    for ind, the_line in enumerate(gtf):  # Parse the entries
        # Parse the attributes column, column 9, of each line
        # noinspection PyTypeChecker
        entry = dict(  # Get the result as dictionary
            [[k, l.replace("\"", "")] for k, l in  # 4) Remove double quote from value
             [j.strip().split(" ", 1) for j in  # 3) For each attribute split key from value
              the_line[8].split(';')  # 1) Split first with ';'
              if j]])  # 2) See if element contains anything (last character is ';', to remove the artifact of it)
        entry.update(dict(zip(columns, the_line[:8])))  # Append dictionary with remaining information (Info in columns)
        # noinspection PyTypeChecker
        if entry["feature"] == "CDS":
            output.append([entry[ft] if ft in entry else np.nan for ft in filter_columns])
        if verbose and (ind % 100 == 0 or ind == len(gtf) - 1):  # Show the progress if 'verbose' is True
            progress_bar(ind, len(gtf) - 1)
    return pd.DataFrame(output, columns=filter_columns)


def download_gtf_refseq(temp_repo_dir, data_url=None):
    """
    Function is to download or find the file in temp_repo_dir.
    :param temp_repo_dir: String. Directory to download or find the file
    :return: String. Path of gtf or gtf.gz
    """
    if not data_url:
        data_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/" \
                   "GCF_000001405.38_GRCh38.p12_genomic.gtf.gz"
    gtf_file = os.path.basename(data_url)  # File name in the FTP server
    # Paths to download the database or look for the database if it already exists
    gtf_gz_path = os.path.join(temp_repo_dir, gtf_file)  # Compressed .gz format
    gtf_path = os.path.splitext(gtf_gz_path)[0]  # Not compressed
    if os.access(gtf_path, os.R_OK) or os.path.isfile(gtf_path):  # Check if the file exist
        return gtf_path  # Return the path if it exist
    else:  # Download otherwise
        subprocess.run((f"cd {temp_repo_dir}; "
                        f"curl -L -R -O {data_url}; "
                        f"gzip -d {gtf_file}"), shell=True)
        assert os.path.isfile(gtf_path) and os.access(gtf_path, os.R_OK)
        return gtf_path  # Return the compressed file


# End
