library(biomaRt)

args = commandArgs(trailingOnly=TRUE)
biomart_version = args[1]
output_path = args[2]

mart = useEnsembl(biomart = "ensembl",  dataset = "hsapiens_gene_ensembl", version = biomart_version)

genes = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "external_synonym"), mart = mart)
genes[genes == ""] = NA

write.table(genes, output_path, row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
