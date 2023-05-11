library(tidyverse)
source("dataset.R")
#source("common.R")


main = function() {
    #annotate_whole()
    annotate_glu()
    annotate_gaba()
    #annotate_subglu()
}

annotate_whole = function(output_dir = "analysis/whole") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster_name,
        1, "SubCTX_1",
        2, "Glu_1",
        3, "GLU_2",
        4, "GLU_3",
        5, "GABA_1",
        6, "GLU_4",
        7, "GLU_5",
        8, "SubCTX_2",
        9, "GLU_6",
        10, "GLU_7",

    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

annotate_glu = function(output_dir = "analysis/glu") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster_name,
        1, "L4_5_IT",
        2, "L6_IT",
        3, "L5_IT",
        4, "L2_3_IT",
        5, "L6_CT",
        6, "L5_ET",
        7, "NP",

    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

annotate_gaba = function(output_dir = "analysis/gaba") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster, ~subclass,
        1, "gaba_1","gaba_1",
        2, "gaba_2","gaba_2",
        3, "gaba_3","gaba_3",
        4, "gaba_4","gaba_4",
        5, "gaba_5","gaba_5",
     
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}
#99 clusters total, including 4 low quality clusters
annotate_subglu = function(output_dir = "analysis/subglu") {
    subclass_annotation = read_csv("analysis/glu/cluster_annotation.csv")
    cluster_annotation = tribble(
        ~cluster_id, ~cluster, ~subclass,~notes,
        "BMA_BLA_1","BMA_BLA 1","BMA_BLA"," ",
        "BMA_BLA_2","BMA_BLA 2","BMA_BLA"," ",
        "BMA_BLA_3","BMA_BLA 3","BMA_BLA"," ",
        "BMA_BLA_4","BMA_BLA 4","BMA_BLA"," ",

        "CA1_1","CA1 1","CA1"," ",
        "CA1_2","CA1 2","CA1"," ",
        "CA1_3","CA1 3","CA1"," ",
        "CA1_4","CA1 4","CA1"," ",
        "CA1_5","CA1 5","CA1"," ",
        "CA1_6","CA1 6","CA1"," ",
        
        "Ca2_CA3_FC_1","CA2 1","CA2"," ",
        "Ca2_CA3_FC_2","CA3 1","CA3"," ",
        "Ca2_CA3_FC_3","FC 1","FC"," ",
        "Ca2_CA3_FC_4","CA3 2","CA3"," ",
        "Ca2_CA3_FC_5","Low Qual 1","Low Qual"," ",
        "Ca2_CA3_FC_6","CA3 3","CA3"," ",
        
        "Car3_1","Car3 1","Car3"," ",
        "Car3_2","Car3 2","Car3"," ",
        "Car3_3","Car3 3","Car3"," ",
        "Car3_4","Car3 4","Car3"," ",
        "Car3_5","Car3 5","Car3"," ",
        "Car3_6","Car3 6","Car3"," ",
        "Car3_7","Car3 7","Car3"," ",
        "Car3_8","Car3 8","Car3"," ",
        "Car3_9","Car3 9","Car3"," ",
        
        "DG_1","DG 1","DG"," ",
        "DG_2","DG 2","DG"," ",
        "DG_3","DG 3","DG"," ",
        "DG_4","DG 4","DG"," ",
        "DG_5","DG 5","DG"," ",
        
        "L2_3_IT_1","L2_3_IT 1","L2_3_IT"," ",
        "L2_3_IT_2","L2_3_IT 2","L2_3_IT"," ",
        "L2_3_IT_3","L2_3_IT 3","L2_3_IT"," ",
        "L2_3_IT_4","L2_3_IT 4","L2_3_IT"," ",
        "L2_3_IT_5","L2_3_IT 5","L2_3_IT"," ",
        "L2_3_IT_6","L2_3_IT 6","L2_3_IT"," ",
        "L2_3_IT_7","L2_3_IT 7","L2_3_IT"," ",

        "L4_5_IT_1","L4_5_IT 1","L4_5_IT"," ",
        "L4_5_IT_2","L4_5_IT 2","L4_5_IT"," ",
        "L4_5_IT_3","L4_5_IT 3","L4_5_IT"," ",
        "L4_5_IT_4","L4_5_IT 4","L4_5_IT"," ",
        "L4_5_IT_5","L4_5_IT 5","L4_5_IT"," ",
        "L4_5_IT_6","L4_5_IT 6","L4_5_IT"," ",
        "L4_5_IT_7","L4_5_IT 7","L4_5_IT"," ",
        
        "L5_ET_1","L5_ET 1","L5_ET"," ",
        "L5_ET_2","L5_ET 2","L5_ET"," ",
        "L5_ET_3","L5_ET 3","L5_ET"," ",
        "L5_ET_4","L5_ET 4","L5_ET"," ",
        "L5_ET_5","L5_ET 5","L5_ET"," ",
        "L5_ET_6","L5_ET 6","L5_ET"," ",
        "L5_ET_7","L5_ET 7","L5_ET"," ",
        "L5_ET_8","L5_ET 8","L5_ET"," ",
        
        "L5_IT_1","L5_IT 1","L5_IT"," ",
        "L5_IT_2","L5_IT 2","L5_IT"," ",
        "L5_IT_3","BLA 1","BLA"," ",
        "L5_IT_4","L5_IT 4","L5_IT"," ",
        "L5_IT_5","L5_IT 5","L5_IT"," ",
        "L5_IT_6","L5_IT 6","L5_IT"," ",
        
        "L6_CT_1","L6_CT 1","L6_CT"," ",
        "L6_CT_2","L6_CT 2","L6_CT"," ",
        "L6_CT_3","L6_CT 3","L6_CT"," ",
        "L6_CT_4","L6_CT 4","L6_CT"," ",
        "L6_CT_5","L6_CT 5","L6_CT"," ",
        "L6_CT_6","L6_CT 6","L6_CT"," ",
        "L6_CT_7","L6_CT 7","L6_CT"," ",
        "L6_CT_8","L6_CT 8","L6_CT"," ",
        
        "L6_IT_1","L6_IT 1","L6_IT"," ",
        "L6_IT_2","L6_IT 2","L6_IT"," ",
        "L6_IT_3","L6_IT 3","L6_IT"," ",
        "L6_IT_4","L6_IT 4","L6_IT"," ",
        "L6_IT_5","L6_IT 5","L6_IT"," ",
        
        "L6b_1","Pir 1","Pir"," ",
        "L6b_2","Pir 2","Pir"," ",
        "L6b_3","L6b 1","L6b"," ",
        "L6b_4","Pir 3","Pir"," ",
        "L6b_5","L6b 2","L6b"," ",
        "L6b_6","L6b 3","L6b"," ",
        "L6b_7","Low_Qual 2","Low_qual"," ",
        
        "Low_Qual_1","Low_Qual 2","Low_Qual"," ",
        "Low_Qual_2","Low_Qual 3","Low_Qual"," ",
        "Low_Qual_3","Low_Qual 4","Low_Qual"," ",
        "Low_Qual_4","Low_Qual 5","Low_Qual"," ",
        "Low_Qual_5","Low_Qual 6","Low_Qual"," ",
        "Low_Qual_6","Low_Qual 7","Low_Qual"," ",
        
        "NP_1","NP 1","NP"," ",
        "NP_2","NP 2","NP"," ",
        "NP_3","NP 3","NP"," ",
        "NP_4","NP 4","NP"," ",
        "NP_5","NP 5","NP"," ",
        "NP_6","NP 6","NP"," ",
        
        "Pir_DL_1","Pir 1","Pir"," ",
        "Pir_DL_2","Pir 2","Pir"," ",
        "Pir_DL_3","Pir 3","Pir"," ",
        "Pir_DL_4","Pir 4","Pir"," ",
        "Pir_DL_5","Pir 5","Pir"," ",
        "Pir_DL_6","Pir 6","Pir"," ",
        
        "Pir_UL_1","Pir 8","Pir"," ",
        "Pir_UL_2","Pir 9","Pir"," ",
        "Pir_UL_3","Pir 10","Pir"," ",
        "Pir_UL_4","Pir 11","Pir"," ",
        "Pir_UL_5","Pir 12","Pir"," ",
        "Pir_UL_6","Pir 13","Pir"," ",
        "Pir_UL_7","Low Qual 8","Low_Qual"," ",
        "Pir_UL_8","Low Qual 9","Low_Qual"," ",
        
        "RSP_UL_1","RSP 1","RSP"," ",
        "RSP_UL_2","RSP 2","RSP"," ",
        "RSP_UL_3","RSP 3","RSP"," ",
        "RSP_UL_4","RSP 4","RSP"," ",
        "RSP_UL_5","RSP 5","RSP"," ",
        "RSP_UL_6","RSP 6","RSP"," ",
        
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

if (sys.nframe() == 0) {
    main()
}

