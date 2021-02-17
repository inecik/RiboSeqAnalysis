"""

"""

from infrastructure.main import *
import time


temp_repo_dir = "/home/raf_pc/Kemal/Temp/human"
data_repo_dir = "/home/raf_pc/Kemal/Data/sam_bukau"
script_path_infrastructure = "/home/raf_pc/Kemal/RiboSeqAnalysis/infrastructure/"

spt = [os.path.join(data_repo_dir, i) for i in ["Sixtymers_TT1.sam", "Sixtymers_TT2.sam"]]
sps = [os.path.join(data_repo_dir, i) for i in ["Sixtymers_Rep1.sam", "Sixtymers_Rep2.sam", "Sixtymers_NoPK.sam"]]

erb_serb = ["/home/raf_pc/Kemal/Data/sam_bukau/SeRP/EBP1/Rep1/IP/IP1.sam",
            "/home/raf_pc/Kemal/Data/sam_bukau/SeRP/EBP1/Rep2/IP/IP2.sam"]
erb_total = ["/home/raf_pc/Kemal/Data/sam_bukau/SeRP/EBP1/Rep1/Total/Total1.sam",
             "/home/raf_pc/Kemal/Data/sam_bukau/SeRP/EBP1/Rep2/Total/Total2.sam"]

nac_serb = ["/home/raf_pc/Kemal/Data/sam_bukau/SeRP/NAC/Rep1/IP/IP1.sam",
            "/home/raf_pc/Kemal/Data/sam_bukau/SeRP/NAC/Rep2/IP/IP2.sam"]
nac_total = ["/home/raf_pc/Kemal/Data/sam_bukau/SeRP/NAC/Rep1/Total/Total1.sam",
             "/home/raf_pc/Kemal/Data/sam_bukau/SeRP/NAC/Rep2/Total/Total2.sam"]

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
                  verbose=True,
                  serb=[["selectiveErb", erb_total, erb_serb], ["selectiveNac", nac_total, nac_serb]],
                  sixtymers=[spt, sps],
                  #coco=[coco_m, coco_d]
                  )


