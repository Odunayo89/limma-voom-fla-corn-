library(edgeR)
library(limma)

counts <- read.csv("fla-corn.csv", row.names = 1, check.names = FALSE)

metadata <- data.frame(
  sample = c("LIC_1","LIC_3","LIC_31","LIC_35","LIC_37",
             "LIF_22","LIF_24","LIF_26","LIF_28","LIF_30"),
  group  = c("LIC","LIC","LIC","LIC","LIC",
             "LIF","LIF","LIF","LIF","LIF"),
  stringsAsFactors = FALSE
)

# 3) Align metadata to counts columns (VERY IMPORTANT)
metadata <- metadata[match(colnames(counts), metadata$sample), ]
stopifnot(all(metadata$sample == colnames(counts)))

# 4) Define group and set baseline
group <- factor(metadata$group)
group <- relevel(group, ref = "LIC")

# 5) limma-voom pipeline (your original steps)
dge <- DGEList(counts = counts, group = group)
dge <- calcNormFactors(dge)

design <- model.matrix(~ group)
v <- voom(dge, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

res_limma <- topTable(fit, coef = 2, number = Inf)
head(res_limma)

# Optional: save results
write.csv(res_limma, "limma_voom_results.csv")

library(ggplot2)
library(dplyr)
library(pheatmap)

# 4) Volcano plot
# ---------------------------
vol <- res_limma %>%
  mutate(gene = rownames(res_limma),
         neglog10FDR = -log10(adj.P.Val),
         sig = adj.P.Val < 0.05 & abs(logFC) > 1)

top_labels <- vol %>% arrange(adj.P.Val) %>% head(10)

p_vol <- ggplot(vol, aes(x = logFC, y = neglog10FDR, color = sig)) +
  geom_point(alpha = 0.7) +
  geom_text(data = top_labels, aes(label = gene),
            hjust = -0.2, vjust = 0.5, size = 3) +
  labs(title = "Volcano plot (limma-voom): LIF vs LIC",
       x = "log2 Fold Change", y = "-log10(FDR)")

print(p_vol)

ggsave("limma_volcano.png", plot = p_vol, width = 10, height = 8, dpi = 300)

# ---------------------------
# 5) MA plot (limma style)
# ---------------------------
ma <- vol %>% mutate(A = AveExpr, M = logFC)

p_ma <- ggplot(ma, aes(x = A, y = M, color = sig)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0) +
  labs(title = "MA plot (limma-voom): LIF vs LIC",
       x = "Average Expression (AveExpr)", y = "log2 Fold Change")

print(p_ma)

ggsave("limma_MA.png", plot = p_ma, width = 10, height = 8, dpi = 300)

# ---------------------------
# 6) Heatmap (top 30 genes by FDR)
# ---------------------------
top30 <- rownames(res_limma[order(res_limma$adj.P.Val), ])[1:30]

mat <- v$E[top30, ]             # voom logCPM
mat_z <- t(scale(t(mat)))       # z-score per gene
mat_z[is.na(mat_z)] <- 0        # safety

ann <- data.frame(group = group)
rownames(ann) <- colnames(mat_z)

pheatmap(mat_z,
         annotation_col = ann,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Top 30 genes heatmap (limma-voom)")

# Save heatmap
png("limma_heatmap_top30.png", width = 1100, height = 900)
pheatmap(mat_z,
         annotation_col = ann,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Top 30 genes heatmap (limma-voom)")
dev.off()


