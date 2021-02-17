"""

"""

from infrastructure.main import *
import time


temp_repo_dir = "/home/raf_pc/Kemal/Temp/human"
data_repo_dir = "/home/raf_pc/Kemal/Data/sam_bukau"
script_path_infrastructure = "/home/raf_pc/Kemal/RiboSeqAnalysis/infrastructure/"

spt = [os.path.join(data_repo_dir, i) for i in ["Sixtymers_TT1.sam", "Sixtymers_TT2.sam"]]
sps = [os.path.join(data_repo_dir, i) for i in ["Sixtymers_Rep1.sam", "Sixtymers_Rep2.sam", "Sixtymers_NoPK.sam"]]

erb_serb = [f"/home/raf_pc/Kemal/Data/sam_bukau/SeRP/EBP1/Rep{i+1}/IP/IP{i+1}.sam" for i in range(2)]
erb_total = [f"/home/raf_pc/Kemal/Data/sam_bukau/SeRP/EBP1/Rep{i+1}/Total/Total{i+1}.sam" for i in range(2)]
nac_serb = [f"/home/raf_pc/Kemal/Data/sam_bukau/SeRP/NAC/Rep{i+1}/IP/IP{i+1}.sam" for i in range(2)]
nac_total = [f"/home/raf_pc/Kemal/Data/sam_bukau/SeRP/NAC/Rep{i+1}/Total/Total{i+1}.sam" for i in range(2)]

coco_d = [os.path.join(data_repo_dir, i) for i in ["Coco_Dis1.sam", "Coco_Dis2.sam"]]
coco_m = [os.path.join(data_repo_dir, i) for i in ["Coco_Mono1.sam", "Coco_Mono2.sam"]]

exclude_genes = ["ENSG00000160789"]

I = Infrastructre(temp_repo_dir,
                  exclude_genes=exclude_genes,
                  riboseq_assign_at=-15,
                  riboseq_assign_to="best_transcript",
                  ensembl_release=102,
                  organism="homo_sapiens",
                  include_gene3d=True,
                  verbose=True)

# Should be first to calculate, since multiprocessing is quite memory inefficient.
I.riboseq_coco = RiboSeqCoco(I.temp_repo_dir, coco_m, coco_d, "cocoassembly",
                             I.riboseq_assign_at, I.riboseq_assign_to,
                             I.protein_genome, I.gene_info,
                             exclude_genes=I.exclude_genes, verbose=I.verbose)

I.riboseq_sixtymers = RiboSeqSixtymers(I.temp_repo_dir, sps, spt, "sixtymers",
                                       I.riboseq_assign_at, I.riboseq_assign_to,
                                       I.protein_genome, I.gene_info,
                                       exclude_genes=I.exclude_genes, verbose=I.verbose)

I.riboseq_selErb = RiboSeqSelective(I.temp_repo_dir, erb_total, erb_serb, "selectiveErb",
                                    I.riboseq_assign_at, I.riboseq_assign_to,
                                    I.protein_genome, I.gene_info,
                                    exclude_genes=I.exclude_genes, verbose=I.verbose)

I.riboseq_selNac = RiboSeqSelective(I.temp_repo_dir, nac_total, nac_serb, "selectiveErb",
                                    I.riboseq_assign_at, I.riboseq_assign_to,
                                    I.protein_genome, I.gene_info,
                                    exclude_genes=I.exclude_genes, verbose=I.verbose)

