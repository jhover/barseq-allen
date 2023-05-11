library(tidyverse)
library(hdf5r)
library(SingleCellExperiment)


convert_20230302 = function() {
  sce = convert_v7_filtneurons("filt_bc_neurons.mat")
  saveRDS(sce, "barseq_20230302.rds")    
}

convert_20230110 = function() {
  sce = convert_v7_filtneurons("filt_neurons-Feb7_newclust.mat")
  saveRDS(sce, "barseq_20230110.rds")    
}


convert_20220822 = function() {
  sce = convert_v7_filtneurons("filt_neurons-bc-RG-removedG-v73.mat")
  saveRDS(sce, "barseq_20220822.rds")    
}


convert_20220621 = function() {
  sce = convert_v7_filtneurons("filt_neurons-bc-RG-fixed-v73.mat")
  saveRDS(sce, "barseq_20220621.rds")    
}


convert_20220127 = function() {
  sce = convert_v7("alldata20220228-bc.mat")
  saveRDS(sce, "barseq_20220228-bc.rds")    
}

convert_220118 = function() {
    sce = convert_v7("alldata20220118-M1plus.mat")
    saveRDS(sce, "barseq_220118.rds")    
}

convert_220102 = function() {
    sce = convert_v7("alldata20220102-M1plus.mat")
    saveRDS(sce, "barseq_220102.rds")    
}

convert_210627 = function() {
    sce = convert_v7("XC199Ldata-20210630.mat")
    saveRDS(sce, "barseq_210630.rds")    
    add_aligned_pos()
    add_ccf("barseq_210630.rds")
}

convert_210525 = function() {
    sce = convert_v7("XC119Ldata.mat")
    saveRDS(sce, "barseq_210525.rds")    
}

convert_210424 = function() {
    all_data = readMat("alldata20210424.mat", fixNames = FALSE)
    sce = convert_v6(all_data)
    saveRDS(sce, "barseq_210424.rds")
}

convert_v7_filtneurons = function(mat_file) {
  h5_file = H5File$new(mat_file)
  #h5_file$ls()
  
  n_genes = length(readDataSet(h5_file[["filt_neurons/expmat/jc"]]))-1
  gene_refs = readDataSet(h5_file[["filt_neurons/genes"]])
  gene_ids = gene_refs$dereference()
  genes = lapply(gene_ids, readDataSet)
  genes = lapply(genes, intToUtf8)
  genes = unlist(genes[1:n_genes])
  
  sample_name = as.character(readDataSet(h5_file[["filt_neurons/id"]]))
  position = readDataSet(h5_file[["filt_neurons/pos"]])
  depth = readDataSet(h5_file[["filt_neurons/depth"]])
  angle = readDataSet(h5_file[["filt_neurons/angle"]])
  slice = readDataSet(h5_file[["filt_neurons/slice"]])
  fov_position = readDataSet(h5_file[["filt_neurons/pos40x"]])
  
  expr = Matrix::sparseMatrix(
    i=readDataSet(h5_file[["filt_neurons/expmat/ir"]])+1,
    p=readDataSet(h5_file[["filt_neurons/expmat/jc"]]),
    x=as.vector(readDataSet(h5_file[["filt_neurons/expmat/data"]])),
    dims=c(length(sample_name), length(genes)),
    dimnames = list(sample_name, genes)
  )
  expr = t(expr)
  
  metadata = data.frame(
    slice = as.vector(slice),
    pos_x = position[,1],
    pos_y = position[,2],
    fov_x = fov_position[,1],
    fov_y = fov_position[,2],
    angle = as.vector(angle),
    depth_x = depth[,1],
    depth_y = depth[,2]
  )
  sce = SingleCellExperiment(assays = list(counts = expr), colData = metadata)
  return(sce)
}
convert_v7 = function(mat_file) {
    h5_file = H5File$new(mat_file)
    #h5_file$ls()

    n_genes = length(readDataSet(h5_file[["neurons/expmat/jc"]]))-1
    gene_refs = readDataSet(h5_file[["neurons/genes"]])
    gene_ids = gene_refs$dereference()
    genes = lapply(gene_ids, readDataSet)
    genes = lapply(genes, intToUtf8)
    genes = unlist(genes[1:n_genes])

    sample_name = as.character(readDataSet(h5_file[["neurons/id"]]))
    position = readDataSet(h5_file[["neurons/pos"]])
    depth = readDataSet(h5_file[["neurons/depth"]])
    angle = readDataSet(h5_file[["neurons/angle"]])
    slice = readDataSet(h5_file[["neurons/slice"]])
    fov_position = readDataSet(h5_file[["neurons/pos40x"]])

    expr = Matrix::sparseMatrix(
        i=readDataSet(h5_file[["neurons/expmat/ir"]])+1,
        p=readDataSet(h5_file[["neurons/expmat/jc"]]),
        x=as.vector(readDataSet(h5_file[["neurons/expmat/data"]])),
        dims=c(length(sample_name), length(genes)),
        dimnames = list(sample_name, genes)
    )
    expr = t(expr)
    
    metadata = data.frame(
        slice = as.vector(slice),
        pos_x = position[,1],
        pos_y = position[,2],
        fov_x = fov_position[,1],
        fov_y = fov_position[,2],
        angle = as.vector(angle),
        depth_x = depth[,1],
        depth_y = depth[,2]
    )
    sce = SingleCellExperiment(assays = list(counts = expr), colData = metadata)
    return(sce)
}

convert_v6 = function(all_data) {
    library(R.matlab)
    
    #rownames(all_data$neurons)
    expr = t(all_data$neurons[,,]$expmat)
    sample_name = as.character(all_data$neurons[,,]$id)
    position = all_data$neurons[,,]$pos
    depth = all_data$neurons[,,]$depth
    angle = all_data$neurons[,,]$angle
    slice = all_data$neurons[,,]$slice
    fov_position = all_data$neurons[,,]$pos40x
    genes = as.character(unlist(all_data$neurons[,,]$genes))
    
    rownames(expr) = genes[1:nrow(expr)]
    colnames(expr) = as.character(sample_name)
    metadata = data.frame(
        slice = as.vector(slice),
        pos_x = position[,1],
        pos_y = position[,2],
        fov_x = fov_position[,1],
        fov_y = fov_position[,2],
        angle = as.vector(angle),
        depth_x = depth[,1],
        depth_y = depth[,2]
    )
    sce = SingleCellExperiment(assays = list(counts = expr), colData = metadata)
    return(sce)
}

convert_v6_filtneurons = function(all_data) {
  library(R.matlab)
  
  #rownames(all_data$filt_neurons)
  expr = t(all_data$filt_neurons[,,]$expmat)
  sample_name = as.character(all_data$filt_neurons[,,]$id)
  position = all_data$filt_neurons[,,]$pos
  depth = all_data$filt_neurons[,,]$depth
  angle = all_data$filt_neurons[,,]$angle
  slice = all_data$filt_neurons[,,]$slice
  fov_position = all_data$filt_neurons[,,]$pos40x
  genes = as.character(unlist(all_data$filt_neurons[,,]$genes))
  
  rownames(expr) = genes[1:nrow(expr)]
  colnames(expr) = as.character(sample_name)
  metadata = data.frame(
    slice = as.vector(slice),
    pos_x = position[,1],
    pos_y = position[,2],
    fov_x = fov_position[,1],
    fov_y = fov_position[,2],
    angle = as.vector(angle),
    depth_x = depth[,1],
    depth_y = depth[,2]
  )
  sce = SingleCellExperiment(assays = list(counts = expr), colData = metadata)
  return(sce)
}
add_aligned_pos = function() {
    barseq = readRDS("barseq_210630.rds")
    pos_matrix = cbind(barseq$pos_x, barseq$pos_y,1)
    rot_matrices = convert_slicetform()
    slice = barseq$slice
    new_pos = matrix(NA, nrow(pos_matrix), ncol(pos_matrix))
    for (i in 1:max(slice)) {
        is_slice = slice == i
        new_pos[is_slice,] = pos_matrix[is_slice,] %*% rot_matrices[[i]]
    }
    any(is.na(new_pos))
    barseq$aligned_x = new_pos[,1]
    barseq$aligned_y = new_pos[,2]
    saveRDS(barseq, "barseq_210630.rds")
}

convert_slicetform = function() {
    all_data = readMat("slicetform1.mat") 
    result = lapply(all_data$slicetform1, "[[",1)
}

add_ccf = function(filename) {
    barseq = readRDS(filename)
    ccf = read_csv("neurons_registered.csv") %>%
        select(sample_id=id, ccf_name=name, ccf_red=red, ccf_green=green, ccf_blue=blue,
               ccf_x = x_CCF, ccf_y = y_CCF, ccf_z = z_CCF)
    ccf$sample_id = as.character(ccf$sample_id)
    barseq$sample_id = colnames(barseq)
    new_coldata = left_join(as.data.frame(colData(barseq)), ccf, by="sample_id")
    colData(barseq) = DataFrame(new_coldata)
    colnames(barseq) = barseq$sample_id
    saveRDS(barseq, filename)
}

if (sys.nframe() == 0) {
  convert_20230302()
}

