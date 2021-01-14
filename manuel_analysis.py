"""

"""


from infrastructure.main import *


temp_repo_dir = "/home/raf_pc/Kemal/Temp"
script_path_infrastructure = "/home/raf_pc/Kemal/RiboSeqAnalysis/infrastructure/"

spt = [os.path.join(temp_repo_dir, i) for i in ["TT1.sam", "TT2.sam"]]
sps = [os.path.join(temp_repo_dir, i) for i in ["Rep1.sam", "Rep2.sam", "NoPK.sam"]]
exclude_genes = np.array(["ENSG00000160789"])


erb_serb = ["/home/raf_pc/Kemal/Data/SeRP/EBP1/Rep1/IP/IP1.sam", "/home/raf_pc/Kemal/Data/SeRP/EBP1/Rep2/IP/IP2.sam"]
erb_total = ["/home/raf_pc/Kemal/Data/SeRP/EBP1/Rep1/Total/Total1.sam", "/home/raf_pc/Kemal/Data/SeRP/EBP1/Rep2/Total/Total2.sam"]
nac_serb = ["/home/raf_pc/Kemal/Data/SeRP/NAC/Rep1/IP/IP1.sam", "/home/raf_pc/Kemal/Data/SeRP/NAC/Rep2/IP/IP2.sam"]
nac_total = ["/home/raf_pc/Kemal/Data/SeRP/NAC/Rep1/Total/Total1.sam", "/home/raf_pc/Kemal/Data/SeRP/NAC/Rep2/Total/Total2.sam"]

I = Infrastructre(temp_repo_dir, exclude_genes=exclude_genes, include_gene3d=True,
                  serb=[["erb", erb_total, erb_serb], ["nac", nac_total, nac_serb]],
                  sixtymers=[spt, sps])




