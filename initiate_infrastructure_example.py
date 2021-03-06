"""

"""


#from infrastructure.main import *


#temp_repo_dir = "/home/raf_pc/Kemal/Temp/human"
#data_repo_dir = "/home/raf_pc/Kemal/Data/sam_bukau"
#script_path_infrastructure = "/home/raf_pc/Kemal/RiboSeqAnalysis/infrastructure/"

# Geçici
temp_repo_dir = "/home/kai/KEMALINECIK/from_raf_computer/Kemal/Temp/human"
data_repo_dir = "/home/kai/KEMALINECIK/from_raf_computer/Kemal/Data/sam_bukau"
script_path_infrastructure = "/home/kai/KEMALINECIK/from_raf_computer/Kemal/RiboSeqAnalysis/infrastructure"
import sys
sys.path.insert(0, '/home/kai/KEMALINECIK/from_raf_computer/Kemal/RiboSeqAnalysis')
from infrastructure.main import *
# Geçici end

spt = [os.path.join(data_repo_dir, i) for i in ["Sixtymers_TT1.sam", "Sixtymers_TT2.sam"]]
sps = [os.path.join(data_repo_dir, i) for i in ["Sixtymers_Rep1.sam", "Sixtymers_Rep2.sam", "Sixtymers_NoPK.sam"]]

erb_serb = [os.path.join(data_repo_dir, f"SeRP/EBP1/Rep{i+1}/IP/IP{i+1}.sam") for i in range(2)]
erb_total = [os.path.join(data_repo_dir, f"SeRP/EBP1/Rep{i+1}/Total/Total{i+1}.sam") for i in range(2)]
nac_serb = [os.path.join(data_repo_dir, f"SeRP/NAC/Rep{i+1}/IP/IP{i+1}.sam") for i in range(2)]
nac_total = [os.path.join(data_repo_dir, f"SeRP/NAC/Rep{i+1}/Total/Total{i+1}.sam") for i in range(2)]

coco_d = [os.path.join(data_repo_dir, i) for i in ["Coco_Dis1.sam", "Coco_Dis2.sam"]]
coco_m = [os.path.join(data_repo_dir, i) for i in ["Coco_Mono1.sam", "Coco_Mono2.sam"]]

exclude_genes = ["ENSG00000160789"]

I = Infrastructre(temp_repo_dir,
                  exclude_genes=exclude_genes,
                  ensembl_release=102,
                  organism="homo_sapiens",
                  #include_gene3d=True,
                  verbose=True)

# Should be first to calculate, since multiprocessing is quite memory inefficient.
# I.riboseq_coco = RiboSeqCoco(I.temp_repo_dir, coco_m, coco_d, "cocoassembly",
#                              riboseq_assign_to="best_transcript", riboseq_assign_at=-15,
#                              protein_genome_instance=I.protein_genome,
#                              gene_info_dictionary=I.gene_info,
#                              exclude_genes=I.exclude_genes, verbose=I.verbose)
#
I.riboseq_sixtymers = RiboSeqSixtymers(I.temp_repo_dir, spt, sps, "sixtymers",
                                        riboseq_assign_to="best_transcript", riboseq_assign_at="auto",
                                        protein_genome_instance=I.protein_genome,
                                        gene_info_dictionary=I.gene_info,
                                        exclude_genes=I.exclude_genes, verbose=I.verbose)

# I.riboseq_selErb = RiboSeqSelective(I.temp_repo_dir, erb_total, erb_serb, "selectiveErb",
#                                     riboseq_assign_to="best_transcript", riboseq_assign_at=-15,
#                                     protein_genome_instance=I.protein_genome,
#                                     gene_info_dictionary=I.gene_info,
#                                     exclude_genes=I.exclude_genes, verbose=I.verbose)
#
# I.riboseq_selNac = RiboSeqSelective(I.temp_repo_dir, nac_total, nac_serb, "selectiveErb",
#                                     riboseq_assign_to="best_transcript", riboseq_assign_at=-15,
#                                     protein_genome_instance=I.protein_genome,
#                                     gene_info_dictionary=I.gene_info,
#                                     exclude_genes=I.exclude_genes, verbose=I.verbose)

