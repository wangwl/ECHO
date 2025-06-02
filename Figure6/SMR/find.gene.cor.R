library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
genes = fread("/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/ref/gencode.v19.genes.hg38.bed", header = FALSE)
colnames(genes) = c("chr", "start", "end", "gene", "score", "strand")
#epi = fread("/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/DICE/filtered/B_CELL_NAIVE.epi", header = FALSE) %>% select(V2)
file = "/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi"
# epi = fread(args[1], header = FALSE) %>% select(V2)
epi = fread(file, header = FALSE) %>% select(V2)
colnames(epi) = c("gene")

output = merge(epi, genes, all.x = T) %>% select(gene, chr, start, strand, score)
output$probe = output$gene
output$chr = gsub("chr", "", output$chr)

output$chr[is.na(output$chr)] <- "Y"
output$chr[grepl("^Un", output$chr)] <- "Y"
output$score[is.na(output$score)] <- 0
output$start[is.na(output$start)] <- 100
output$strand[is.na(output$strand)] <- "+"

setcolorder(output, c("chr", "probe", "score", "start", "gene", "strand"))

write.table(output, file = "/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi.gene.cor.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# write.table(output, file = args[2], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# write.table(output, file = "/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/DICE/filtered/B_CELL_NAIVE.gene.cor.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
