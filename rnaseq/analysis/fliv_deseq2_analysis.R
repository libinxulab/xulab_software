# fliv_deseq2_analysis.R
# --- Dylan H. Ross
# --- 2021/08/29
# --- Description:
# --- --- Performs differential axpression analysis of Vanessa's female mouse liver RNAseq data from 
# --- --- BAC-C16 and control mice using DESeq2. 
# --- Requires:
# --- --- C16_CTRL_fliv_counts.csv (count data)
# --- --- individual_genes (a directory to store plots of individual gene counts)
# --- Produces:
# --- --- plots of mean counts vs. LFC with and without LFC shrinkage applied
# --- --- plots of individual gene counts (stored in individual_genes directory)
# --- --- LFC and statistical significance data (.csv format, filtered and unfiltered by p-value)
# --- --- variance stabilizing transformed gene count data (.csv format)


# load DESeq2 library, ignore startup messages
suppressMessages(library("DESeq2"))
# used for LFC shrinkage (Zhu, Ibrahim, and Love 2018)
library("apeglm")
# enable parallelism
library("BiocParallel")
register(MulticoreParam(4))


# counts_df_from_csv
# --- loads the raw count data from "BC16_CTRL_fliv_counts_raw.csv"
# --- puts the data into a data frame
# --- removes any row with sum of 0 counts across samples
# --- returns the data frame
counts_df_from_csv <- function() {
    df <- data.frame(read.csv("BC16_CTRL_fliv_counts_raw.csv", row.names = "Geneid"))
    return(df[apply(df[,-1], 1, function(x) !all(x == 0)),])  # return df with zero sum rows removed
}


# init_deseq_data
# --- initializes a DESeqDataSet instance from counts dataframe
# --- returns initialized DESeqDataSet
init_deseq_data <- function(counts_df) {
    # counts matrix
    cts <- as.matrix(counts_df) 
    # sample names with corresponding experimental condition (BC16 vs. CTRL)
    coldata <- data.frame(row.names = colnames(counts_df), condition = c(rep('BC16', 10), rep('CTRL', 8)))
    coldata$condition <- as.factor(coldata$condition) 
    # experimental design is to compare by condition (BC16 vs. CTRL)
    return(DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition))
} 


# plotMA_
# --- wrapper around plotMA function
# --- saves plot as .png image 
plotMA_ <- function(deseq_res, alpha, ylim, png) {
    png(png, 
        width = 3, 
        height = 3, 
        units = "in", 
        res = 600, 
        pointsize = 4)
    plotMA(deseq_res, alpha = alpha, ylim = ylim)
    dev.off()
}


# plot_indiv_gene_counts
plot_indiv_gene_counts <- function(deseq_data, gene) {
    png(paste0("individual_genes/", gene, "_counts.png"),
        width = 1.25, 
        height = 1.75, 
        units = "in", 
        res = 600, 
        pointsize = 4)
    d <- plotCounts(deseq_data, gene = gene, intgroups = "condition", returnData = TRUE)
    x <- split(d$count, d$condition)
    boxplot(x, xlab = "condition", ylab = "counts", medlwd = 1.5, main = gene)
    points(jitter(rep(1, length(x[[1]])), factor = 5), 
           x[[1]], pch = 1, lwd = 1., cex = 1.5)
    points(jitter(rep(2,length(x[[2]])), factor = 5), 
           x[[2]], pch = 1, lwd = 1., cex = 1.5)
    dev.off()
}


# main
# --- main execution sequence
main <- function() {

    # load raw cound data from csv and initialize DESeq data instance
    deseq_data <- init_deseq_data(counts_df_from_csv())

    # perform differential expression analysis
    deseq_data <- DESeq(deseq_data)
    print(head(deseq_data))
    
    # gather results
    deseq_res <- results(deseq_data)
    summary(deseq_res)
    print(mcols(deseq_res)$description)

    # plot mean count vs. LFC (no shrinkage)
    plotMA_(deseq_res, alpha = 0.05, ylim = c(-4, 4), png = "MA_plot_no_shrinkage.png")

    # export results to .csv (no shrinkage)
    write.csv(as.data.frame(deseq_res), 
              file = "fliv_deseq2_results.csv")
    write.csv(subset(as.data.frame(deseq_res), padj < 0.05), 
              file = "fliv_deseq2_results_padj<0.05.csv")

    # apply LFC shrinkage (Zhu, Ibrahim, and Love 2018)
    deseq_res_norm <- lfcShrink(deseq_data, coef = 2, type = "apeglm")
    summary(deseq_res_norm)
    print(mcols(deseq_res_norm)$description)

    # plot some individual gene counts
    print("plotting individual_genes (padj < 0.05)...")
    for (gene in rownames(subset(as.data.frame(deseq_res_norm), padj < 0.05))) {
        plot_indiv_gene_counts(deseq_data, gene = gene)
    }
    print("... done")

    # plot mean count vs. LFC (shrinkage applied)
    plotMA_(deseq_res_norm, alpha = 0.05, ylim = c(-4, 4), png = "MA_plot_shrinkage.png")

    # export results to .csv (shrinkage applied)
    write.csv(as.data.frame(deseq_res_norm), 
              file = "fliv_deseq2_results_shrinkage.csv")
    write.csv(subset(as.data.frame(deseq_res_norm), padj < 0.05), 
              file = "fliv_deseq2_results_shrinkage_padj<0.05.csv")

    # export vst transformed counts
    write.csv(assay(vst(deseq_data, blind = FALSE, nsub = 5000)), 
              file = "fliv_deseq2_counts_vst.csv")

}

main()
