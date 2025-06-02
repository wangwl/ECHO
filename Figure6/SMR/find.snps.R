library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

snps = fread("/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/smr/cis-meQTL/cis-meQTL_snps.txt", header = FALSE)
colnames(snps) = c("chr", "pos", "snp", "ref", "alt")

# esi = fread("/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/DICE/filtered/B_CELL_NAIVE.esi", header = FALSE)
# file = "/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi"
# esi = fread(file, header = FALSE)
esi = fread(args[1], header = FALSE)
colnames(esi) = c("chr", "snp", "0", "pos", "ref", "alt", "maf")

output = merge(snps, esi, by = "snp", sort=FALSE, all.y = TRUE) %>% select(snp, chr.x, pos.x, ref.x, alt.x)
output$chr.x = gsub("chr", "", output$chr.x)
output$dist = 0

# output[is.na(output$chr.x)]
output$chr.x[is.na(output$chr.x)] <- "Y"
output$pos.x[is.na(output$pos.x)] <- 100
output$ref.x[is.na(output$ref.x)] <- "T"
output$alt.x[is.na(output$alt.x)] <- "A"

setcolorder(output, c("chr.x", "snp", "dist", "pos.x", "ref.x", "alt.x"))

# write.table(output, file = "/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/DICE/filtered/B_CELL_NAIVE.snp.cor.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(output, file = "/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi.snp.cor.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# write.table(output, file = args[2], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
