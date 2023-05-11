
library(tidyverse)
library(SingleCellExperiment)

BARSEQ_DIR = ""


load_barseq = function(filename="barseq_20220920.rds", slice = NULL, normalization_factor = 10, min_genes = 5, min_counts = 20) {
    #sce = readRDS(file.path(BARSEQ_DIR, filename))
    sce = readRDS(filename)
    rownames(sce) = convert_gene_names(rownames(sce))
    sce = sce[, colSums(counts(sce)) >= min_counts & colSums(counts(sce)>0) >= min_genes]
    colnames(sce) = paste0(sce$slice, "_", colnames(sce))
    cpm(sce) = convert_to_cpm(counts(sce), normalization_factor)
    if (!is.null(slice)) {
        sce = sce[, sce$slice %in% slice]
    }
    return(sce)
}

convert_gene_names = function(genes) {
    to_convert = c("Tafa1" = "Fam19a1", "Tafa2" = "Fam19a2", "Ccn2"="Ctgf")
    needs_conversion = genes %in% names(to_convert)
    genes[needs_conversion] = to_convert[genes[needs_conversion]]
    return(genes)
}

convert_to_cpm = function(M, total_counts = 1000000) {
    normalization_factor = Matrix::colSums(M) / total_counts
    if (is(M, "dgCMatrix")) {
        M@x = M@x / rep.int(normalization_factor, diff(M@p))
        return(M)
    } else {
        return(scale(M, center = FALSE, scale = normalization_factor))
    }
}

read_panel = function(panel_name) {
    read_csv(file.path("panel", paste0(panel_name, ".csv"))) %>%
        select(cell_type, gene)
}

subclass_markers = function() {
    marker_panel = read_panel("final_panel") %>%
        select(cell_type, gene) %>%
        filter(!(cell_type %in% c("IT", "IT RSP", "IT TPE-ENT", "IT/PT", "PT"))) %>%
        mutate(group = "all")
    return(marker_panel)
}

previous_panel = function() {
    my_genes = unique(read_panel("test_panel")$gene)
    my_genes = c(my_genes, "Slc17a7", "Gad1", "Slc30a3")
    return(my_genes)
}

read_clusters = function() {
    read_csv("analysis/cluster.csv")
}
