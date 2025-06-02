library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# input = "/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/smr/output/eQTLGen"
# smr <- fread(paste0(input, "/B-Mem_eQTLGen.smr"))

file = args[1]
smr <- fread(file)

# filter out the HLA region 
MHC_start = 28477797
MHC_end = 33448354 

smr = smr %>% filter(!(Outco_Chr == 6 & Outco_bp < MHC_end & Outco_bp > MHC_start))

# calculate the bonferroni threshold 
bonferroni = (0.05 / nrow(smr))

# signif = smr %>% filter(p_SMR < bonferroni) %>% select(Expo_ID, Outco_ID, p_SMR)
signif = smr %>% filter(p_SMR < bonferroni)

# output = "/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/smr/output/eQTLGen"
# write.table(signif, file = paste0(output, "/B-Mem_eQTLGen.signif.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

write.table(signif, file = args[2], quote = FALSE, sep = "\t", row.names = FALSE)

num_dmrs = signif$Expo_ID %>% unique() %>% length()
print(paste0("number of significant DMRs: ", num_dmrs))
percent = num_dmrs/ (smr$Expo_ID %>% unique() %>% length())
print(paste0("percent of DMRs with at least one significant association: ", percent))

print(paste0("bonferroni threshold:", bonferroni))

