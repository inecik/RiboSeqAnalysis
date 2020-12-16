library(biomaRt)

args = commandArgs(trailingOnly=TRUE)
biomart_version = args[1]
output_path = args[2]

mart = useEnsembl(biomart = "ensembl",  dataset = "hsapiens_gene_ensembl", version = biomart_version)

genes = getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id",
                           "uniprotswissprot", "uniprotsptrembl", 
                           "external_gene_name", "external_synonym", 
                           "transcript_biotype",
                           "transcript_tsl", "transcript_gencode_basic", "transcript_appris", "transcript_mane_select", "external_transcript_name",
                           "chromosome_name", "start_position", "end_position", "strand"
                           ), mart = mart)
genes[genes == ""] = NA

print(output_path)
write.table(genes, output_path, row.names=FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
