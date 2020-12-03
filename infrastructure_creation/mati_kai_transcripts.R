# Kemal Inecik

args = commandArgs(trailingOnly=TRUE)
r_data_path = args[1]
output_path = args[2]

load(r_data_path)

base1=sigmoidfits[["HEK"]][["1"]][["base"]][["transcript"]]
ssig1=sigmoidfits[["HEK"]][["1"]][["ssig"]][["transcript"]]
dsig1=sigmoidfits[["HEK"]][["1"]][["dsig"]][["transcript"]]

base2=sigmoidfits[["HEK"]][["2"]][["base"]][["transcript"]]
ssig2=sigmoidfits[["HEK"]][["2"]][["ssig"]][["transcript"]]
dsig2=sigmoidfits[["HEK"]][["2"]][["dsig"]][["transcript"]]

transcripts = c(base1, ssig1, dsig1, base2, ssig2, dsig2)

write.table(transcripts, output_path, row.names=FALSE, col.names=FALSE, quote = FALSE, sep = "\n")

