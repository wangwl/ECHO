library(ggplot2)
library(data.table)
library(gridExtra)
library(dplyr)

dir = "/Users/rosanwang/Documents/research/rotations/ecker/echo/inputs"

window_start = 65.6 * 10^6
# window_start = min(min(eQTL$pos), min(mQTL$pos))
window_end = 66 * 10^6
# window_end = max(max(eQTL$pos), max(mQTL$pos))

EFEMP2_eQTL = fread(paste0(dir, "/eQTLGen_EFEMP2.tsv")) %>% 
  select(SNPChr, SNPPos_hg38, `-log10p`) %>% 
  filter(SNPPos_hg38 > window_start, SNPPos_hg38 < window_end)
colnames(EFEMP2_eQTL) = c("chr", "pos", "-log10p")
EFEMP2_eQTL$label = "eQTL"
EFEMP2_eQTL$color = "eQTL"

eQTL = EFEMP2_eQTL

mQTL = fread(paste0(dir, "/meQTL_dmr223620.tsv")) %>% 
  select(var_chr, var_from, `-log10p`) %>% 
  filter(var_from > window_start, var_from < window_end)
colnames(mQTL) = c("chr", "pos", "-log10p")
mQTL$label = "mQTL"
mQTL$color = "mQTL"

gwas = fread(paste0(dir, "/GCST90044763.cc_no_beta_tsv")) %>% 
  select(chrom, position, pvalues) %>% 
  filter(chrom == "chr11", position > window_start, position < window_end)
gwas$pvalues = log10(gwas$pvalues) * -1 
colnames(gwas) = c("chr", "pos", "-log10p")
gwas$color = "gwas"
gwas$label = "gwas"

# filter for snps that are tested in all three 
common_values <- Reduce(intersect, list(gwas$pos, mQTL$pos, eQTL$pos))

QTL = rbind(gwas, eQTL, mQTL)
# QTL = filter(QTL, QTL$pos %in% common_values)
# dmr223620: 65791784
# p value = 1.424159e-08

# rs10791824 - chr11:65791795
QTL$color[QTL$pos == 65791795] <- "top_gwas"

QTL$label = factor(QTL$label, levels=c('gwas','eQTL','mQTL'))

smr_bonferroni = -1 * log10(1.79670821494523e-08)

plot = ggplot() + geom_point(data = QTL, aes(x = pos, y = `-log10p`, color = color)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_manual(values = c("top_gwas" = "red", "eQTL" = "maroon", "mQTL" = "darkblue", "gwas" = "grey"))

hline_df <- data.frame(label = "gwas", yint = smr_bonferroni)

# Add hline only to group A panel
plot = plot + geom_hline(data = hline_df, aes(yintercept = yint), color = "red", linetype = "dashed")

# add gene information
# genes = fread(paste0(dir, "/glist-hg38"))
genes = fread(paste0(dir, "/EFEMP2.loci.txt"))
colnames(genes) = c("chr", "start", "end","strand", "geneID", "geneType", "geneName")
genes$geneType = sub('.*"(.*)".*', '\\1', genes$geneType)
genes$geneName = sub('.*"(.*)".*', '\\1', genes$geneName)

genes = genes %>% filter(chr == "chr11", start > window_start, 
                         end < window_end, geneType == "protein_coding")
genes$y = seq(0, by = 0.75, length.out = nrow(genes))
genes$start_copy = genes$start
genes$start <- ifelse(genes$strand == "-", genes$start, genes$end)
genes$end <- ifelse(genes$strand == "-", genes$end, genes$start_copy) 
genes$label = "annot"

# add smr probe pval
smr_pval = 5.047654e-10
new_row <- data.frame(chr = 11, pos = 65791784, `-log10p` = log10(smr_pval) * -1, label = "gwas", color = "mQTL", check.names = FALSE)
plot = plot + geom_point(data = new_row, aes(x = pos, y = `-log10p`), shape = 18, size = 4, color = "blue")


plot <- plot + 
  geom_segment(data = genes, aes(x = end, y = y, xend = start, yend = y),
               size = 1, arrow = arrow(length = unit(0.1, "inches"))) +
  geom_text(data = genes, aes(x = start, y = y, label = geneName),
            hjust = 'outside', nudge_x = -1) + 
  facet_grid(rows = vars(label), scales = "free_y")


plot



# ggsave("EFEMP2_loci.pdf", plot = plot)















# plot(eQTL$SNPPos_hg38, eQTL$`-log10p`, pch = 16)
# plot(mQTL$var_from, mQTL$`-log10p`, pch = 16)
# plot(eQTL$SNPPos_hg38, eQTL$`-log10p`, ylim=c(y.min,y.max), xlim=c(x.min,x.max),
# cex.axis=axis,xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
