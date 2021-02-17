library(biomaRt)

args = commandArgs(trailingOnly=TRUE)
biomart_version = args[1]
dataset_name = args[2]
output_path = args[3]

mart = useEnsembl(biomart = "ensembl",  dataset = dataset_name, version = biomart_version)

genes = getBM(attributes=c("ensembl_gene_id", "uniprotswissprot", "uniprotsptrembl"), mart = mart)
genes[genes == ""] = NA

write.table(genes, output_path, row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
