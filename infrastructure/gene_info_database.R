library(biomaRt)

args = commandArgs(trailingOnly=TRUE)
biomart_version = args[1]
dataset_name = args[2]
output_path = args[3]

mart = useEnsembl(biomart = "ensembl",  dataset = dataset_name, version = biomart_version)

att = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "transcript_biotype",
                     "transcript_tsl", "transcript_gencode_basic", "transcript_appris", "transcript_mane_select",
                     "external_transcript_name", "chromosome_name", "start_position", "end_position", "strand")
if (dataset_name != "hsapiens_gene_ensembl") {
  att = att[att != "transcript_mane_select"]
}

genes = getBM(attributes=att, mart = mart)
genes[genes == ""] = NA
write.table(genes, output_path, row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
