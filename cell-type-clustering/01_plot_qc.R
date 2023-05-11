
library(tidyverse)
source("dataset.R")


main = function() {
    # current QC criteria
    min_genes = 5
    min_counts = 20

    output_dir = file.path("figs", "QC")
    dir.create(output_dir,recursive = TRUE, showWarnings = FALSE)
    barseq = load_barseq(min_genes = 0, min_counts = 0)
    qc = data.frame(
        sample = colnames(barseq),
        n_genes = colSums(counts(barseq) > 0),
        n_counts = colSums(counts(barseq))
    )

    qc %>%
        ggplot(aes(x = n_genes)) +
        geom_bar() +
        theme_classic(base_size = 20) +
        geom_vline(xintercept = 5, linetype="dashed")
    ggsave(file.path(output_dir, "n_genes.pdf"))
    
    qc %>%
        ggplot(aes(x = n_counts)) +
        geom_histogram(bins = 50) +
        theme_classic(base_size = 20) +
        geom_vline(xintercept = 20, linetype="dashed")
    ggsave(file.path(output_dir, "n_counts.pdf"))
    
    f = round(100*mean(qc$n_genes >= min_genes & qc$n_counts >= min_counts))
    n = round(sum(qc$n_genes >= min_genes & qc$n_counts >= min_counts)/1e6,2)
    qc %>%
        filter(n_genes > 0) %>%
        ggplot(aes(x = n_genes, y = n_counts)) +
        ggrastr::geom_jitter_rast(alpha=0.01) +
        geom_density_2d() +
        theme_classic(base_size=20) +
        scale_x_log10() +
        scale_y_log10() +
        geom_hline(yintercept = min_counts, linetype = "dashed") +
        geom_vline(xintercept = min_genes, linetype = "dashed") +
        ggtitle(paste0(n, "M (", f, "%) cells pass QC criteria"))
    ggsave(file.path(output_dir, "counts_vs_genes.pdf"))
}

if (sys.nframe() == 0) {
    main()
}

