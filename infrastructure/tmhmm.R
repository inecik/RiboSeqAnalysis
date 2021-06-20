library(biomaRt)

args = commandArgs(trailingOnly=TRUE)
biomart_version = args[1]
dataset_name = args[2]
output_path = args[3]

mart = useEnsembl(biomart = "ensembl",  dataset = dataset_name, version = biomart_version)

tmhmm = getBM(attributes=c("ensembl_peptide_id", "tmhmm_start", "tmhmm_end"), mart = mart)
tmhmm[tmhmm == ""] = NA
tmhmm = na.omit(tmhmm)

write.table(tmhmm, output_path, row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
