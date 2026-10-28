args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript PlotHostPhageSummary.R known.tsv positive.tsv output_dir min_species_cells min_host_cells main_host_contribution")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggbeeswarm)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

known <- read.delim(args[1], colClasses = "character")
positive <- read.delim(args[2], colClasses = "character")
result_dir <- args[3]
min_species_cells <- as.integer(args[4])
min_host_cells <- as.integer(args[5])
main_host_contribution <- as.numeric(args[6])
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot_obj, filename, width, height) {
  ggsave(file.path(result_dir, filename), plot = plot_obj, width = width, height = height)
}

species_total <- aggregate(Cell ~ SGB, unique(known), length)
names(species_total)[2] <- "AllCellNum"
species_positive <- aggregate(Cell ~ SGB, unique(positive[c("Cell", "SGB")]), length)
names(species_positive)[2] <- "PhageCellNum"
species <- merge(species_total, species_positive, by = "SGB", all.x = TRUE)
species$PhageCellNum[is.na(species$PhageCellNum)] <- 0L
species$PhageCellProp <- species$PhageCellNum / species$AllCellNum
species <- species[species$AllCellNum > min_species_cells, ]
if (nrow(species) < 2) stop("Fewer than two SGBs pass min_species_cells")

cell_phages <- aggregate(MGroup ~ Cell, positive, function(x) length(unique(x)))
frequency <- as.data.frame(table(cell_phages$MGroup))
names(frequency) <- c("PhageNum", "CellNum")
p_step2 <- ggplot(frequency, aes(PhageNum, CellNum)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = CellNum), vjust = -0.5, size = 4) +
  labs(title = "Phage Types in Cell", x = "Number of Phages Types per Cell", y = "CellNum") +
  scale_y_continuous(limits = c(0, 9000), breaks = seq(0, 8000, by = 1000)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
save_plot(p_step2, "phage_types_per_cell.pdf", 7, 5)

set.seed(123)
p_step4 <- ggplot(species, aes(x = "", y = PhageCellProp)) +
  geom_violin(width = 0.8, fill = "#3498DB", alpha = 0.8, trim = FALSE) +
  geom_quasirandom(aes(size = AllCellNum), method = "pseudorandom", width = 0.4,
                   alpha = 0.7, shape = 21, fill = "white", color = "black") +
  scale_fill_viridis_d(option = "D", alpha = 0.3) +
  scale_size_continuous(range = c(1, 8), name = "Cell Count") +
  stat_summary(fun = median, geom = "point", shape = 18, size = 5,
               color = "#F08080", alpha = 0.8) +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)),
               vjust = -1, color = "black", size = 4) +
  theme_minimal() +
  labs(title = "Violin: The PhageCell Proportion in each SGB") +
  guides(size = guide_legend(override.aes = list(alpha = 1)))
save_plot(p_step4, "phage_positive_species_proportion.pdf", 6, 6)

p_step5 <- ggplot(species, aes(AllCellNum, PhageCellNum)) +
  geom_point(alpha = 0.7, size = 4, shape = 21, color = "white", fill = "#3498DB") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E74C3C", fill = "#F1948A") +
  ggpubr::stat_cor(method = "pearson",
                   aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
                   label.x.npc = "left", label.y.npc = "top", size = 5) +
  ggpubr::stat_regline_equation(label.x.npc = "left", label.y.npc = 0.9, size = 5) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  labs(title = "Relationship Between Host SGB AllCellNum and PhageCellNum Variables",
       subtitle = "With linear regression and correlation analysis",
       x = "AllCellNum", y = "PhageCellNum") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 12))
save_plot(p_step5, "phage_positive_cell_correlation.pdf", 7, 6)

host_pairs <- positive %>%
  group_by(MGroup, SGB) %>%
  summarize(Count = n(), .groups = "drop")
host_counts <- host_pairs %>%
  pivot_wider(names_from = MGroup, values_from = Count,
              values_fill = list(Count = 0)) %>%
  as.data.frame()
rownames(host_counts) <- host_counts$SGB
host_counts <- as.matrix(host_counts[, -1, drop = FALSE])
host_counts <- host_counts[rowSums(host_counts) > min_species_cells,
                           colSums(host_counts) > min_species_cells, drop = FALSE]
host_counts <- host_counts[rownames(host_counts) != "no", , drop = FALSE]
host_counts <- host_counts[, colSums(host_counts) != 0, drop = FALSE]
if (nrow(host_counts) == 0 || ncol(host_counts) == 0) stop("No phage-host matrix passes the cell thresholds")
host_share <- sweep(host_counts, 2, colSums(host_counts), "/")

# Match the reference's cumulative squared contribution rule for main hosts.
main_pairs <- do.call(rbind, lapply(colnames(host_counts), function(mgroup) {
  cells <- sort(host_counts[, mgroup][host_counts[, mgroup] > 0], decreasing = TRUE)
  stop_at <- which(cumsum(cells^2) / sum(cells^2) >= main_host_contribution)[1]
  data.frame(MGroup = mgroup, SGB = names(cells)[seq_len(stop_at)],
             CellNum = as.integer(cells[seq_len(stop_at)]))
}))
major_pairs <- main_pairs[main_pairs$CellNum >= min_host_cells, ]
if (nrow(major_pairs) == 0) stop("No major phage-host pair passes min_host_cells")
major_mgroups <- unique(major_pairs$MGroup)
bar_data <- data.frame(MGroup = major_mgroups,
                       CellNum = colSums(host_counts[, major_mgroups, drop = FALSE]))
bar_data <- bar_data[order(bar_data$CellNum), ]
bar_data$MGroup <- factor(bar_data$MGroup, levels = bar_data$MGroup)

host_plot_matrix <- host_share[unique(main_pairs$SGB),
                               as.character(bar_data$MGroup), drop = FALSE]
host_tree <- pheatmap::pheatmap(
  host_plot_matrix, cluster_cols = FALSE, border_color = "gray50",
  color = colorRampPalette(c("#9db3ef", "white", "firebrick"))(100),
  breaks = seq(min(host_plot_matrix), max(host_plot_matrix), length.out = 100),
  silent = TRUE
)$tree_row
heatmap_data <- t(host_plot_matrix)
marks <- matrix("", nrow = nrow(heatmap_data), ncol = ncol(heatmap_data),
                dimnames = dimnames(heatmap_data))
shown_pairs <- main_pairs[main_pairs$MGroup %in% rownames(marks), ]
marks[cbind(match(shown_pairs$MGroup, rownames(marks)),
            match(shown_pairs$SGB, colnames(marks)))] <- "*"
host_number <- as.data.frame(table(major_pairs$MGroup))
names(host_number) <- c("MGroup", "HostNum")
host_number <- host_number[match(as.character(bar_data$MGroup), host_number$MGroup), ]
host_annotations <- ComplexHeatmap::rowAnnotation(
  `Cell Number` = ComplexHeatmap::anno_barplot(
    log2(bar_data$CellNum + 1), gp = grid::gpar(fill = "orange", col = NA),
    which = "row", border = FALSE, bar_width = 0.6, width = grid::unit(2.5, "cm"),
    axis_param = list(at = c(0, 5, 10), labels_rot = 0)
  ),
  `Host Number` = ComplexHeatmap::anno_barplot(
    host_number$HostNum, baseline = 0, ylim = c(0, 12), which = "row",
    gp = grid::gpar(fill = "#5eab81", col = NA), border = FALSE,
    bar_width = 0.6, width = grid::unit(2.5, "cm"),
    axis_param = list(at = seq(0, 12, by = 2), labels_rot = 0)
  ),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  gap = grid::unit(2, "mm")
)
host_figure <- ComplexHeatmap::Heatmap(
  heatmap_data, name = "Cell fraction",
  heatmap_legend_param = list(direction = "horizontal", legend_width = grid::unit(4, "cm")),
  col = circlize::colorRamp2(
    c(min(heatmap_data), mean(range(heatmap_data)), max(heatmap_data)),
    c("#9db3ef", "white", "firebrick")
  ),
  rect_gp = grid::gpar(col = "gray50", lwd = 0.3),
  cluster_rows = FALSE, cluster_columns = as.dendrogram(host_tree),
  column_dend_reorder = FALSE, show_column_dend = TRUE,
  column_dend_height = grid::unit(1.5, "cm"),
  column_dend_gp = grid::gpar(col = "black", lwd = 0.5),
  right_annotation = host_annotations,
  row_names_side = "left", row_names_gp = grid::gpar(fontsize = 7),
  column_names_side = "top", column_names_rot = 90,
  column_names_gp = grid::gpar(fontsize = 7),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (marks[i, j] == "*") grid::grid.text("*", x, y, gp = grid::gpar(fontsize = 8))
  }
)
pdf(file.path(result_dir, "phage_host_range_summary.pdf"), width = 16, height = 18)
ComplexHeatmap::draw(host_figure, heatmap_legend_side = "bottom")
dev.off()
unlink(file.path(result_dir, c("phage_host_cell_counts.pdf", "phage_host_breadth.pdf")))
