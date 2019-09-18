args = commandArgs(trailingOnly=TRUE)

options(stringsAsFactors = F)

VarScan <- read.table(args[1], skip = 22, sep = "\t")

diff_w <- VarScan$V5 * VarScan$V6
VarScan.filter <- VarScan[VarScan$V5 == 0 | VarScan$V6 == 0 | diff_w < 0,]

Gwas_info <- read.table("Data/Variants/Gwas_v1.0_hg19.bed", sep = "\t", quote = "\"")[,c(4,7,8)]

VarScan.filter <- merge(VarScan.filter, Gwas_info, by.x = 'V2', by.y = 'V4')
colnames(VarScan.filter) <- c("var_id","ac_motif", "var_class","var_coord","best_w", "worst_w", "w_diff", "best_pval", "worst_pval", "pval_ratio", "best_variant", 
  "worst_variant", "best_offset", "worst_offset", "min_offset_diff", "best_strand", "worst_strand", "str_change",
  "best_seq", "worst_seq", "minor_allele_freq", "Gwas_disease", "Gwas_gene")

write.table(VarScan.filter, file = args[2], quote = F, sep = "\t", row.names = F)
