# 基础统计
sgb_name <- "$1"

# 读取数据
cells <- read.table("StrainCells.txt", header=TRUE, sep="\t")
snp_data <- read.table(paste0(sgb_name, "_SNPpd.txt"), header=TRUE, sep="\t", check.names=FALSE)

cat("=== 基础统计结果 ===\n")
cat("SGB名称:", sgb_name, "\n")
cat("总细胞数:", nrow(cells), "\n")
cat("菌株数:", length(unique(cells$Cluster)), "\n")
cat("SNP位点数:", nrow(snp_data), "\n")

# 菌株分布
cat("\n菌株分布:\n")
strain_table <- table(cells$Cluster)
for(strain in names(strain_table)) {
    cat(sprintf("  菌株 %s: %d 个细胞\n", strain, strain_table[strain]))
}

# 保存结果
stats <- data.frame(
    项目 = c("SGB名称", "总细胞数", "菌株数", "SNP位点数"),
    数值 = c(sgb_name, nrow(cells), length(unique(cells$Cluster)), nrow(snp_data))
)
write.table(stats, "basic_statistics.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("\n统计结果已保存到: basic_statistics.tsv\n")
