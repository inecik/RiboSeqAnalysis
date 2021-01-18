"""

"""

from infrastructure.main import *
import time


temp_repo_dir = "/home/raf_pc/Kemal/Temp"
data_repo_dir = "/home/raf_pc/Kemal/Data"
script_path_infrastructure = "/home/raf_pc/Kemal/RiboSeqAnalysis/infrastructure/"

spt = [os.path.join(data_repo_dir, "Sixtymers", i) for i in ["60mer_TT1.sam", "60mer_TT2.sam"]]
sps = [os.path.join(data_repo_dir, "Sixtymers", i) for i in ["60mer_Rep1.sam", "60mer_Rep2.sam", "60mer_NoPK.sam"]]

erb_serb = ["/home/raf_pc/Kemal/Data/SeRP/EBP1/Rep1/IP/IP1.sam", "/home/raf_pc/Kemal/Data/SeRP/EBP1/Rep2/IP/IP2.sam"]
erb_total = ["/home/raf_pc/Kemal/Data/SeRP/EBP1/Rep1/Total/Total1.sam", "/home/raf_pc/Kemal/Data/SeRP/EBP1/Rep2/Total/Total2.sam"]

nac_serb = ["/home/raf_pc/Kemal/Data/SeRP/NAC/Rep1/IP/IP1.sam", "/home/raf_pc/Kemal/Data/SeRP/NAC/Rep2/IP/IP2.sam"]
nac_total = ["/home/raf_pc/Kemal/Data/SeRP/NAC/Rep1/Total/Total1.sam", "/home/raf_pc/Kemal/Data/SeRP/NAC/Rep2/Total/Total2.sam"]

coco_d = [os.path.join(data_repo_dir, "Coco", i) for i in ["Dis1.sam", "Dis2.sam"]]
coco_m = [os.path.join(data_repo_dir, "Coco", i) for i in ["Mono1.sam", "Mono2.sam"]]


exclude_genes = ["ENSG00000160789"]

I = Infrastructre(temp_repo_dir, exclude_genes=exclude_genes,
                  #include_gene3d=True,
                  #serb=[["selectiveErb", erb_total, erb_serb], ["selectiveNac", nac_total, nac_serb]],
                  #sixtymers=[spt, sps],
                  coco=[coco_m, coco_d])


