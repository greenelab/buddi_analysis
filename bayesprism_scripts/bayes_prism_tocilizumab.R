library(BayesPrism)
require(data.table)
require(tidyr)


read_and_format <- function(in_file){
    
    sc_matr = fread(in_file, header=T)
    cell_ids = colnames(sc_matr)[2:ncol(sc_matr)]
    sc_matr = data.frame(sc_matr)
    
    # format the column names
    colnames(sc_matr)[1] = "hgnc"

    # take all value columns
    count_matr = sc_matr[,2:ncol(sc_matr)]
    
    # format the data for BayesPrism
    gene_ids = sc_matr$hgnc
    count_matr = data.frame(count_matr)
    rownames(count_matr) = gene_ids
    colnames(count_matr) = cell_ids

    return(count_matr)
    
}

write_for_getway <- function(sc_ref, bulk_ref, sc_dir, out_name){


    # make the cell reference matrix
    cell_type = colnames(sc_ref)
    cell_names = 1:ncol(sc_ref)
    cell_ref_df = data.frame(cell_id=cell_names)
    cell_ref_df$cell_type = cell_type
    cell_ref_df$cell_subtype = cell_type
    cell_ref_df$tumor_state = rep(0, length(cell_type))

    # make same cell id names
    colnames(sc_ref) = cell_names

    out_file = paste0(sc_dir, out_name, "_gateway_bulk.tsv")
    write.table(bulk_ref, out_file, sep="\t", quote=F, row.names = T, col.names=T)
    out_file = paste0(sc_dir, out_name, "_gateway_bulk.rds")
    saveRDS(as.data.frame(bulk_ref), out_file)

    out_file = paste0(sc_dir, out_name, "_gateway_scRef.tsv")
    write.table(sc_ref, out_file, sep="\t", quote=F, row.names = T, col.names=T)
    out_file = paste0(sc_dir, out_name, "_gateway_scRef.rds")
    saveRDS(as.data.frame(sc_ref), out_file)

    out_file = paste0(sc_dir, out_name, "_gateway_cellRef.csv")
    write.table(cell_ref_df, out_file, sep=",", quote=F, row.names = T, col.names=T)
    out_file = paste0(sc_dir, out_name, "_gateway_cellRef.rds")
    saveRDS(as.data.frame(cell_ref_df), out_file)


}

format_bayes_prism <- function(sc_dir, bulk_dir, file_id_sc, file_id_bulk, out_name){
    sc_matr_file = paste0(sc_dir, file_id_sc, ".tsv")
    sc_ref = read_and_format(sc_matr_file)
    
    bulk_matr_file = paste0(bulk_dir, file_id_bulk, ".tsv")
    bulk_ref = read_and_format(bulk_matr_file)


    write_for_getway(sc_ref, bulk_ref, sc_dir, out_name)

    
}

format_bp_Z <- function(cell_type_name, bp_res_Z){

    grep_str = paste0("[.]", cell_type_name, "$")
    ct_idx = grep(grep_str, colnames(bp_res_Z))
    ct_z = bp_res_Z[,ct_idx]
    colnames(ct_z) = substr(colnames(ct_z), 1, nchar(colnames(ct_z))-(nchar(cell_type_name)+1))
    print(dim(ct_z)[2] == 13999)

    ct_meta = data.frame(cell_type=rep(cell_type_name, dim(ct_z)[1]))
    ct_meta$full_name = row.names(ct_z)
    ct_meta = separate(data = ct_meta, col = full_name, into = c("sample_id", "time"), sep="-")


    return(list(ct_z, ct_meta))


}

read_bayes_prism_results <- function(bp_Z_file, bp_file,  out_dir){

    # get the cell type proportions
    bp_res = readRDS(bp_file)

    bp_res = data.frame(bp_res)
    bp_res$full_name = row.names(bp_res)
    bp_res = separate(data = bp_res, col = full_name, into = c("sample_id", "time"), sep="-")

    outfile = paste0(out_dir, "bp_tocilizumab_prop_01.tsv")
    write.table(bp_res, outfile, sep="\t", quote=F, row.names=F)


    # get the gene expression
    bp_res_Z = readRDS(bp_Z_file)


    res = format_bp_Z("NK", bp_res_Z)
    ct_z = res[[1]]
    ct_meta = res[[2]]


    res = format_bp_Z("T_CD4", bp_res_Z)
    ct_z = rbind(ct_z, res[[1]])
    ct_meta = rbind(ct_meta, res[[2]])

    res = format_bp_Z("T_other", bp_res_Z)
    ct_z = rbind(ct_z, res[[1]])
    ct_meta = rbind(ct_meta, res[[2]])


    res = format_bp_Z("B", bp_res_Z)
    ct_z = rbind(ct_z, res[[1]])
    ct_meta = rbind(ct_meta, res[[2]])


    res = format_bp_Z("T_CD8", bp_res_Z)
    ct_z = rbind(ct_z, res[[1]])
    ct_meta = rbind(ct_meta, res[[2]])


    res = format_bp_Z("Myeloid", bp_res_Z)
    ct_z = rbind(ct_z, res[[1]])
    ct_meta = rbind(ct_meta, res[[2]])


    res = format_bp_Z("Fibroblast", bp_res_Z)
    ct_z = rbind(ct_z, res[[1]])
    ct_meta = rbind(ct_meta, res[[2]])


    res = format_bp_Z("Endothelial", bp_res_Z)
    ct_z = rbind(ct_z, res[[1]])
    ct_meta = rbind(ct_meta, res[[2]])

    outfile = paste0(out_dir, "bp_tocilizumab_Z_01.tsv")
    write.table(ct_z, outfile, sep="\t", quote=F, row.names=F)


    outfile = paste0(out_dir, "bp_tocilizumab_Z_meta_01.tsv")
    write.table(ct_meta, outfile, sep="\t", quote=F, row.names=F)

}


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

sc_dir = args[1]
bulk_dir = args[2]
file_id_sc = args[3]
file_id_bulk = args[4]
bp_Z_file = args[5]
bp_file = as.numeric(args[6])
out_name = as.numeric(args[7])
out_dir = as.numeric(args[8])

print("sc_dir:")
print(sc_dir)

print("bulk_dir:")
print(bulk_dir)

print("file_id_sc:")
print(file_id_sc)

print("file_id_bulk:")
print(file_id_bulk)

print("bp_Z_file:")
print(bp_Z_file)

print("bp_file:")
print(bp_file)

print("out_name:")
print(out_name)

print("out_dir:")
print(out_dir)

read_bayes_prism_results(bp_Z_file, bp_file, out_dir)
format_bayes_prism(sc_dir, bulk_dir, file_id_sc, file_id_bulk, out_name)



