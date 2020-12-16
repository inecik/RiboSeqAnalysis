"""

"""

# DISULFIDE PROBLEMI

import os
import sys
import numpy as np
import joblib
import pysam
from skimage.filters import *
from matplotlib import pyplot as plt
from scipy.ndimage import *
from module_conservation.functions import *

#__file__ = "/home/kai/Ribo-seq-Analysis/module_conservation/uniprot_structure_parser.py" # todo

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")
from module_supplementary.common_functions import progressBarForTerminal
from module_supplementary.common_functions import bcolors as c

feature_index = joblib.load("/Users/kemalinecik/Desktop/genome_base_feature_index_SHORT.joblib")
genome_base = joblib.load("/Users/kemalinecik/Desktop/genome_base_SHORT.joblib")



# Union of features/domains is consistent or not
for ind, i in enumerate(genome_base):
    pass
    footprints = genome_base[i]["footprints"][genome_base[i]['mane_transcript_cds_PA']]
    conservation = genome_base[i]["conservation"][genome_base[i]['mane_transcript_cds_PA']]
    if sum(footprints) > 0:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,2.5))
        smooth_footprints = smooth(footprints, window_len=65, window="hanning")
        smooth_conservation = smooth(conservation, window_len=65*2, window="hanning")
        ax.plot(smooth_footprints/max(smooth_footprints), color="blue")
        ax.plot(footprints/max(footprints), color='gray', alpha=0.5)
        ax.plot(smooth_conservation/max(smooth_conservation), color="red")
        plt.title(f"{i}     Count: {sum(footprints)}")
        ax.figure.tight_layout()
        plt.show()
        #fig.savefig(f"_delete_{i}.pdf")
        plt.close('all')


# domain 0:500
expected = (45/2000) * 500
observed = sum(footprints[0:500])
ratio = (observed-expected)/expected

plt.plot(smooth(genome_base[i]["conservation"], window_len=200));plt.show()

"""
counter = 0
for i in genome_base:
    a=genome_base[i]
    b=a["footprints"][a['mane_transcript_cds_PA']]
    if max(b) > 1:
        d=a["conservation"][a['mane_transcript_cds_PA']]


        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,6))
        ax[0].figsize = (10,2.5)
        ax[1].figsize = (10,2.5)

        counter += 1
        ax[0].plot(b/max(b), color='gray', alpha=0.5)
        ax[0].set_title(f'Footprint: {sum(b)}')
        ax[0].figure.tight_layout()
        ax[1].plot(smooth(d, window_len=100), color='green', alpha=0.5)
        ax[1].set_title('Conservation')
        ax[1].figure.tight_layout()

        fig.suptitle(f"{i}     Length: {len(b)}", y=1.08)
        #plt.show()

        fig.savefig(f"_delete_{i}.pdf")
        plt.close('all')
"""


"""Simple example where a few features are defined "by hand" and are displayed
and exported as PNG, first with a linear view, then with a circular
view.
"""

features = [
    GraphicFeature(
        start=5, end=20, strand=0, color="#ffd700"
    ),
    GraphicFeature(
        start=21,
        end=500,
        strand=-1,
        color="#ffcccc",
        label="Gene 2",
        thickness=20,
    ),
    GraphicFeature(
        start=501, end=701, strand=0, color="#cffccc", label="Gene 2"
    ),
    GraphicFeature(
        start=702, end=900, strand=0, color="#ccccff", label="Gene 3"
    ),
]


from matplotlib import pyplot as plt
from dna_features_viewer import (
    GraphicFeature,
    GraphicRecord)

def graphic_feature_creator(color, ranges, label=None):
    ranges = sorted([sorted(i) for i in ranges])
    return [GraphicFeature(start=i, end=j, color=color, label=label) for i, j in ranges]

def boolean_array_to_ranges(boolean_array):
    ranges = list()
    previous = False
    for ind, element in enumerate(boolean_array):
        if element and not previous:
            ranges.append(ind)
            previous = True
        elif not element and previous:
            ranges.append(ind)
            previous = False
        elif previous and ind == len(boolean_array) - 1:
            ranges.append(ind)
    return [ranges[i:i+2] for i in range(0, len(ranges), 2)]

for gene_id in genome_base:
    annotations = [j for i in ["tracks", "gene3d"] for j in genome_base[gene_id][i]]

    fig, axes = plt.subplots(nrows=len(annotations)+2, ncols=1, figsize=(10, 3 * (len(annotations) + 2)))

    b=genome_base[gene_id]["footprints"][genome_base[gene_id]['mane_transcript_cds_PA']]
    smooth_footprints = smooth(b, window_len=65, window="hanning")
    axes[0].plot(smooth_footprints / max(smooth_footprints), color="blue")
    axes[0].plot(b / max(b), color='gray', alpha=0.5)
    #axes[0].figure.tight_layout()
    if genome_base[gene_id]["strand"] == "-":
        x1, x2 = axes[0].get_xlim()
        axes[0].set_xlim(x2, x1)
    axes[0].spines['right'].set_visible(False)
    axes[0].spines['top'].set_visible(False)

    d = genome_base[gene_id]["conservation"][genome_base[gene_id]['mane_transcript_cds_PA']]
    axes[1].plot(smooth(d, window_len=65), color='red')
    #axes[1].figure.tight_layout()
    if genome_base[gene_id]["strand"] == "-":
        x1, x2 = axes[1].get_xlim()
        axes[1].set_xlim(x2, x1)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)

    i = 2
    for ann_type in ["tracks", "gene3d"]:
        for ann in genome_base[gene_id][ann_type]:
            bool_arr = genome_base[gene_id][ann_type][ann][genome_base[gene_id]['mane_transcript_cds_PA']]
            ranges = boolean_array_to_ranges(bool_arr)
            features = graphic_feature_creator(color=None, ranges=ranges)
            record = GraphicRecord(sequence_length=len(bool_arr), features=features)
            record.plot(ax=axes[i])
            #axes[i].figure.tight_layout()
            axes[i].text(0.5, 0.5, ann, verticalalignment='center', horizontalalignment='center',)

            if genome_base[gene_id]["strand"] == "-":
                x1, x2 = axes[i].get_xlim()
                axes[i].set_xlim(x2, x1)
            i += 1
    fig.savefig(f"_delete_{gene_id}.pdf")
    plt.close("all")


