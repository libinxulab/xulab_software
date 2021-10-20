# two_group_gene_counts_pca.R
# --- Dylan H. Ross
# --- 2021/10/20
# --- --- Computes PCA on gene count data (raw or VST transformed from DESeq2) and plots the projections for the first 
# --- --- three principal components.


# plot_pca_proj
# --- plot the PCA projections for two specified PCs for two groups
# --- takes as input results from prcomp, principal components x and y (e.g. 1 and 2 to plot components 1 and 2), 
# --- and the number of samples in groups A and B
# --- produces PCA_proj_PC<x>_PC<y>.png
plot_pca_proj <- function(pca_res, pc_x, pc_y, n_A, n_B) {
    max_proj <- max(max(pca_res$x[, pc_x]), max(pca_res$x[, pc_y])) + 5
    lim <- c(-max_proj, max_proj)
    pc_x_var <- 100. * summary(pca_res)$importance[2, pc_x]
    pc_y_var <- 100. * summary(pca_res)$importance[2, pc_y]
    png(paste0("PCA_proj_PC", pc_x, "_PC", pc_y, ".png"), 
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
         xlab = paste0("PC", pc_x, sprintf(" (%.1f", pc_x_var), " %)"), 
         ylab = paste0("PC", pc_y, sprintf(" (%.1f", pc_y_var), " %)"), 
         bty = "n")
    points(pca_res$x[1:n_A, pc_x], 
           pca_res$x[1:n_A, pc_y], 
           col = "blue")
    text(pca_res$x[1:n_A, pc_x], 
         pca_res$x[1:n_A, pc_y], 
         rownames(pca_res$x)[1:n_A], 
         pos = 4, offset = 0.25, col = "blue")
    points(pca_res$x[(n_A + 1):(n_B + n_A), pc_x], 
           pca_res$x[(n_A + 1):(n_B + n_A), pc_y], 
           col = "red")
    text(pca_res$x[(n_A + 1):(n_B + n_A), pc_x], 
         pca_res$x[(n_A + 1):(n_B + n_A), pc_y], 
         rownames(pca_res$x)[(n_A + 1):(n_B + n_A)], 
         pos = 4, offset = 0.25, col = "red")
    abline(h = 0, lty = 3, lwd = 0.5)
    abline(v = 0, lty = 3, lwd = 0.5)
    dev.off()
}


# main
# --- main execution sequence
main <- function() {

    # ! SET RUN PARAMETERS HERE !
    # VST counts file
    vst_counts_file <- "CTRL_BC16_counts_raw.csv"
    # number of samples in group A
    n_A <- 4
    # number of samples in group B
    n_B <- 5

    # load count data
    counts_vst <- data.frame(t(read.csv(vst_counts_file, row.names = 1)))

    # compute PCA
    counts_vst_pca <- prcomp(counts_vst)
    print(summary(counts_vst_pca))
    print(counts_vst_pca$x[, 1:3])

    # plot projections
    plot_pca_proj(counts_vst_pca, 1, 2, n_A, n_B)
    plot_pca_proj(counts_vst_pca, 2, 3, n_A, n_B)

}

main()
