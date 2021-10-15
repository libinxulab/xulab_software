# counts_vst_pca.R
# --- Dylan H. Ross
# --- 2021/08/31
# --- --- PCA on VST-transformed gene counts (DESeq2) from  female mouse liver 
# --- --- RNAseq data from BAC-C16 and control mice 


# plot_pca_proj
# --- plot the PCA projections for two specified PCs
plot_pca_proj <- function(pca_res, pc_a, pc_b, padj = FALSE) {
    if (padj) {
        base_name <- "PCA_projections_padj_PC"
        lim <- c(-15, 15)
    } else {
        base_name <- "PCA_projections_PC"
        lim <- c(-65, 65)
    }
    png(paste0(base_name, pc_a, "_PC", pc_b, ".png"), 
        width = 2.5, 
        height = 2.5, 
        units = "in", 
        res = 600, 
        pointsize = 4)
    plot(0, 0, 
         pch = 3,
         lwd = 0.5,
         xlim = lim, 
         ylim = lim,
         xlab = paste0("PC", pc_a), 
         ylab = paste0("PC", pc_b), 
         bty = "n")
    points(pca_res$x[1:10, pc_a], 
           pca_res$x[1:10, pc_b], 
           col = "red")
    text(pca_res$x[1:10, pc_a], 
         pca_res$x[1:10, pc_b], 
         rownames(pca_res$x)[1:10], 
         pos = 4, offset = 0.25, col = "red")
    points(pca_res$x[11:18, pc_a], 
           pca_res$x[11:18, pc_b], 
           col = "blue")
    text(pca_res$x[11:18, pc_a], 
         pca_res$x[11:18, pc_b], 
         rownames(pca_res$x)[11:18], 
         pos = 4, offset = 0.25, col = "blue")
    abline(h = 0, lty = 3, lwd = 0.5)
    abline(v = 0, lty = 3, lwd = 0.5)
    dev.off()
}


# main
# --- main execution sequence
main <- function() {

    # load count data
    counts_vst <- data.frame(t(read.csv("fliv_deseq2_counts_vst.csv", row.names = "Geneid")))

    # compute PCA
    counts_vst_pca <- prcomp(counts_vst)
    print(summary(counts_vst_pca))
    print(counts_vst_pca$x[, 1:3])

    # plot projections
    plot_pca_proj(counts_vst_pca, 1, 2)
    plot_pca_proj(counts_vst_pca, 2, 3)

    # filter counts by padj < 0.05 AND abs(log2fc) >= 1
    sig_genes_present <- c()
    sig_genes <- data.frame(read.csv("fliv_deseq2_results_shrinkage_padj<0.05.csv"))[c("Geneid", "log2FoldChange")]
    for (i in 1 : length(sig_genes$Geneid)) {
        if (is.element(sig_genes$Geneid[i], colnames(counts_vst))) {
            if (abs(sig_genes$log2FoldChange[i]) >= 1) {
                sig_genes_present <- c(sig_genes_present, sig_genes$Geneid[i])
            }
        }
    }
    counts_vst_padj <- counts_vst[, sig_genes_present]

    # compute new PCA with only the significant genes
    counts_vst_padj_pca <- prcomp(counts_vst_padj)
    print(summary(counts_vst_padj_pca))
    print(counts_vst_padj_pca$x[, 1:3])

    # plot projections
    plot_pca_proj(counts_vst_padj_pca, 1, 2, padj = TRUE)
    plot_pca_proj(counts_vst_padj_pca, 2, 3, padj = TRUE)
}

main()
