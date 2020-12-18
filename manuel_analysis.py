"""

"""


from infrastructure.main import *


temp_repo_dir = "/home/raf_pc/Kemal/Data"
spt = [os.path.join(temp_repo_dir, i) for i in ["TT1.sam", "TT2.sam"]]
sps = [os.path.join(temp_repo_dir, i) for i in ["Rep1.sam", "Rep2.sam", "NoPK.sam"]]
exclude_genes = np.array(["ENSG00000160789"])
I = Infrastructre(temp_repo_dir, exclude_genes=exclude_genes, include_gene3d=True)

