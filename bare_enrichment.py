"""

"""

import os
import sys
import numpy as np
import joblib


__file__ = "/home/kai/Ribo-seq-Analysis/module_conservation/uniprot_structure_parser.py" # todo

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")
from module_supplementary.common_functions import progressBarForTerminal
from module_supplementary.common_functions import bcolors as c
from module_conservation.functions import *

# genome_base = joblib.load("/Users/kemalinecik/Desktop/genome_base_SHORT.joblib")
genome_base = joblib.load("/home/kai/KEMALINECIK/out/OUTPUT/protein_structure/genome_base.joblib")


######################## Enrichment
domains = dict()
counter = 0
for gene_id in genome_base:
    mane_bool = genome_base[gene_id]["mane_transcript_cds_PA"]
    footprint_mane = genome_base[gene_id]["footprints"][mane_bool]
    if sum(footprint_mane) > 1:
        counter += 1
        footprint_mane = smooth(footprint_mane, window_len=65)
        for domain_id in genome_base[gene_id]["tracks"]:
            domain_mane = genome_base[gene_id]["tracks"][domain_id][mane_bool]
            # shift 30
            shifted_domain_mane = np.concatenate([domain_mane, np.full(15, False)])
            shifted_footprint_mane = np.concatenate([np.zeros(15), footprint_mane])
            # end
            nucleotide_expected = sum(footprint_mane) / len(shifted_footprint_mane)
            domain_experiment = sum(shifted_footprint_mane[shifted_domain_mane])/(nucleotide_expected * sum(shifted_domain_mane))
            domains[domain_id] = domains.get(domain_id, []) + [domain_experiment]

merged_dom = dict()
for i in domains:
    k='.'.join(i.split('.')[:2])
    merged_dom[k] = merged_dom.get(k, []) + domains[i]


for i in merged_dom:
    lent = len( np.array(merged_dom[i])[~np.isnan(merged_dom[i])] )
    print(f"{i}:\tlength: {lent}\tmedian: {round(np.nanmedian(np.array(merged_dom[i])),2)}")
    if np.nanmedian(np.array(merged_dom[i])) <  0.25 and lent > 100:
        print(i)

seaborn.boxplot(merged_dom["2.90"]); plt.show()
len(merged_dom["2.90"])

for i in domains:
    if np.nanmedian(np.array(domains[i])) :
        print(i)
        print(np.nanmedian(domains[i]))

seaborn.boxplot(domains["LIPID"]); plt.show()


##################### dom Conservation


domains = list()
counter = 0
for gene_id in genome_base:
    mane_bool = genome_base[gene_id]["mane_transcript_cds_PA"]
    if len(genome_base[gene_id]["conservation"][mane_bool]) > 130:
        conservation_mane = smooth(genome_base[gene_id]["conservation"][mane_bool], window_len=30*2)
        footprint_mane = smooth(genome_base[gene_id]["footprints"][mane_bool], window_len=30)
        if not np.isnan(np.nanmean(conservation_mane)) and sum(footprint_mane) > 1:
            expected = np.nanmean(conservation_mane) * sum(footprint_mane)
            observed = sum(conservation_mane*footprint_mane)
            #
            domains.append((observed-expected)/expected)

seaborn.boxplot(domains, showfliers = False); plt.show()


##################### dom Conservation and structure

domains = dict()
counter = 0
for gene_id in genome_base:
    mane_bool = genome_base[gene_id]["mane_transcript_cds_PA"]
    conservation_mane = smooth(genome_base[gene_id]["conservation"][mane_bool], window_len=30 * 2)
    footprint_mane = smooth(genome_base[gene_id]["footprints"][mane_bool], window_len=30)
    if sum(footprint_mane) > 1:
        for domain_id in genome_base[gene_id]["tracks"]:
            domain_mane = genome_base[gene_id]["tracks"][domain_id][mane_bool]
            dom_con = np.nanmean(conservation_mane[domain_mane] * footprint_mane[domain_mane])
            dom_con_exp = np.nanmean(conservation_mane[domain_mane] * sum(footprint_mane) / len(footprint_mane))
            domains[domain_id] = domains.get(domain_id, []) + [dom_con/dom_con_exp]





