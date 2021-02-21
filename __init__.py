xsa = [59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71]
lens = sorted({j for i in I.riboseq_sixtymers.experiment.length_distribution for j in i.keys()})
res = I.riboseq_sixtymers.experiment.assign_certain_lengths(xsa, 40, 0, I.protein_genome, I.gene_info)
calcs = {l: I.riboseq_sixtymers.experiment.periodicity_window("stop", 200, gene_assignments = res[l]["gene"], flanking_assignments = res[l]["flanking"]) for l in res.keys()}

####################################################################################################
k=[]
calcs2 = {lm: np.average(calcs[lm], axis=0, weights=I.riboseq_sixtymers.experiment.total_assigned) for lm in calcs.keys()}
for l in xsa:
    #I.riboseq_sixtymers.experiment.plot_periodicity(calcs[l], "stop", str(l))
    rk = [stats.spearmanr(calcs2[l][:100], calcs2[l][i:100+i])[0] for i in range(19)]
    k.append(rk)
k = np.array(k)

k[k<0.5] = np.nan
plt.gca().imshow(k.T)
plt.gca().set_xticks(range(len(xsa)))
plt.gca().set_yticks(list(range(19)))
# ... and label them with the respective list entries
plt.gca().set_xticklabels(xsa)
plt.gca().set_yticklabels(list(range(19)))
plt.setp(plt.gca().get_xticklabels(), rotation=90, ha="right", va="center",
         rotation_mode="anchor")
plt.show()
####################################################################################################


####################################################################################################
"""
m = []
for l in xsa:
    rk = [stats.spearmanr(calcs[k][0][:80], calcs[l][0][:80])[0] for k in xsa]
    m.append(rk)
plt.gca().imshow(m)
plt.gca().set_xticks(range(len(xsa)))
plt.gca().set_yticks(range(len(xsa)))
# ... and label them with the respective list entries
plt.gca().set_xticklabels(xsa)
plt.gca().set_yticklabels(xsa)
plt.setp(plt.gca().get_xticklabels(), rotation=90, ha="right", va="center",
         rotation_mode="anchor")
plt.show()
"""
####################################################################################################






####################################################################################################
plt.figure( figsize=(14, 4))
for lm in xsa[7:9]:
    k = np.average(calcs[lm], axis=0, weights=I.riboseq_sixtymers.experiment.total_assigned)
    if sum(k) == 0:
        continue
    if any([np.isnan(i) for i in k]):
        continue
    auc, buc = np.nan, np.nan
    k = np.average(calcs[lm], axis=0, weights=I.riboseq_sixtymers.experiment.total_assigned)
    p = k
    v = np.vstack([
        np.add.reduceat(k[0:-3], np.arange(0, len(k) - 3, 3)),
        np.add.reduceat(k[1:-2], np.arange(0, len(k) - 3, 3)),
        np.add.reduceat(k[2:-1], np.arange(0, len(k) - 3, 3))
    ])
    m = v / np.nanmedian(v, axis=1).reshape(v.shape[0], 1)
    a, b = np.where(m > 2.5)
    au, bu = np.unique(a), [sorted(b[a==i])[-1] for i in np.unique(a)]
    auc, buc = au[v[au, bu].argmax()], bu[v[au, bu].argmax()]

    #2, 55
    2 + 55*3 + 0 == 200-33 # stop start
    2 + 55*3 + 2 + 0 == 200-30 # stop end
    o = (200-33 - auc - buc * 3)
    plt.plot(np.arange(len(p)) + o, p,marker=".")
plt.vlines(auc + buc * 3 + o,min(p), max(p), color="red")
plt.vlines(auc + buc * 3 + 2 + o, min(p), max(p), color="red")
plt.show()


    plt.title(f"{lm}")
    plt.show()
    print(f"{lm}: {auc + buc * 3}")
    print(f"Offset: {200-33 - auc - buc * 3}")
    print(f"{lm-(200-33 - auc - buc * 3)}\n")

####################################################################################################


## Average not mean
## always collapse replicates
## filter out significantly small

backup = I.riboseq_sixtymers.experiment.length_distribution

probs = [{k: i[k] / I.riboseq_sixtymers.experiment.total_assigned[ind] for k in i.keys()}
         for ind, i in enumerate(I.riboseq_sixtymers.experiment.length_distribution)]

probs = [{k: i[k] for k in i.keys() if i[k] > 0.025}
         for ind, i in enumerate(probs)]
I.riboseq_sixtymers.experiment.length_distribution = probs
I.riboseq_sixtymers.experiment.plot_footprint_length_distribution()
I.riboseq_sixtymers.experiment.length_distribution = backup

np.array([6583202., 3075599., 2119062.])

"""
- Determine offset from 5'-0_position by above method

"""

# add_riboseq function
I.riboseq_sixtymers_l = RiboSeqSixtymers(I.temp_repo_dir, spt, sps, "sixtymers_l",
                                       riboseq_assign_to="best_transcript", riboseq_assign_at="auto",
                                       protein_genome_instance=I.protein_genome,
                                       gene_info_dictionary=I.gene_info,
                                       exclude_genes=I.exclude_genes, verbose=I.verbose, recalculate=True)
