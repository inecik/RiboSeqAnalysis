"""

"""

from infrastructure.main import *
import time

temp_repo_dir = "/home/raf_pc/Kemal/Temp/mouse"
data_repo_dir = "/home/raf_pc/Kemal/Data/sam_arpat"
script_path_infrastructure = "/home/raf_pc/Kemal/RiboSeqAnalysis/infrastructure/"

disomes = [os.path.join(data_repo_dir, i) for i in ["SRR9715828.sam", "SRR9715826.sam"]]
monosomes = [os.path.join(data_repo_dir, i) for i in ["SRR1930189.sam", "SRR1930188.sam"]]

I = Infrastructre(temp_repo_dir,
                  ensembl_release=102,
                  organism="mus_musculus",
                  include_gene3d=True,
                  verbose=True)

I.riboseq_sixtymers = RiboSeqSixtymers(I.temp_repo_dir, monosomes, disomes, "sixtymers",
                                       selection="best_transcript", assignment=-15,
                                       protein_genome_instance=I.protein_genome,
                                       gene_info_dictionary=I.gene_info,
                                       exclude_genes=I.exclude_genes, verbose=I.verbose,
                                       footprint_len_experiment=list(range(45,71)),  # From paper
                                       footprint_len_translatome=list(range(26,36))  # From paper
                                       )
