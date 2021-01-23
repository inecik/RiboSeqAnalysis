"""
It is to test how much my calculations correlates with the calculations in the paper.
"""

# from infrastructure.main import *

import seaborn as sns

temp_repo_dir = "/home/raf_pc/Kemal/Temp"
data_repo_dir = "/home/raf_pc/Kemal/Data"
script_path_infrastructure = "/home/raf_pc/Kemal/RiboSeqAnalysis/infrastructure/"

coco_d = [os.path.join(data_repo_dir, "Coco", i) for i in ["Dis1.sam", "Dis2.sam"]]
coco_m = [os.path.join(data_repo_dir, "Coco", i) for i in ["Mono1.sam", "Mono2.sam"]]

exclude_genes = ["ENSG00000160789"]

I = Infrastructre(temp_repo_dir, exclude_genes=exclude_genes, coco=[coco_m, coco_d])
paper = pd.read_excel("/home/raf_pc/Kemal/RiboSeqAnalysis/notebooks/coco_correlation_with_paper/abc7151_Table_S1.xlsx", header=1)

comparison = list()
for gene_id in I.gene_list:
    my_onset = I.riboseq_coco.calculate_onset(gene_id)
    for gene_name in I.gene_info[gene_id].gene_names:
        paper_onset = paper[paper["Gene"] == gene_name]["Onset [codons]"]
        if len(paper_onset) != 0:
            comparison.append([gene_id, gene_name, my_onset/3, int(paper_onset)])


comp = pd.DataFrame.from_records(comparison, columns=["gene_id", "gene_name", "mine", "paper"])
splot = sns.relplot("mine", "paper", data=comp)
plt.show()

correlation = comp.corr()
print(correlation)


