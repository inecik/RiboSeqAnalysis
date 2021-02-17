library(biomaRt)

args = commandArgs(trailingOnly=TRUE)
biomart_version = args[1]
dataset_name = args[2]
output_path = args[3]

mart = useEnsembl(biomart = "ensembl",  dataset = dataset_name, version = biomart_version)

gene3d = getBM(attributes=c("ensembl_peptide_id", "gene3d", "gene3d_start", "gene3d_end"), mart = mart)
gene3d[gene3d == ""] = NA
gene3d = na.omit(gene3d)
gene3d = gene3d[, c(1, 3, 4, 2)]

write.table(gene3d, output_path, row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
