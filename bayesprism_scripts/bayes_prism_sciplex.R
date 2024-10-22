library(BayesPrism)
require(data.table)
require(tidyr)


read_and_format <- function(in_file, sc=TRUE){
    
    sc_matr = fread(in_file, header=T)

    # get the cell IDs
    if(sc){
        cell_ids = sc_matr$scpred_CellType
        # take all value columns
        count_matr = sc_matr[,3:ncol(sc_matr)]
    }else{
        cell_ids = paste0("samp", 1:(nrow(sc_matr)))
        # take all value columns
        count_matr = sc_matr[,2:ncol(sc_matr)]
    }
    
    # format the data for BayesPrism
    gene_ids = colnames(count_matr)

    count_matr = t(count_matr)

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

format_bayes_prism <- function(in_dir, out_dir, file_id_train, file_id_test){

    sc_matr_file = paste0(in_dir, file_id_train, "_cybersort_sig.tsv.gz")
    sc_ref = read_and_format(sc_matr_file, sc=T)
    
    bulk_matr_file = paste0(in_dir, file_id_test, "_cybersort_mix.tsv.gz")
    bulk_ref = read_and_format(bulk_matr_file, sc=F)

    write_for_getway(sc_ref, bulk_ref, out_dir, file_id_train)

    
}


run_bp <- function(out_dir, file_id_train, ncores){

    out_file = paste0(out_dir, file_id_train, "_gateway_bulk.rds")
    bk.dat = readRDS(out_file)

    out_file = paste0(out_dir, file_id_train, "_gateway_scRef.rds")
    sc.dat = readRDS(out_file)

    out_file = paste0(out_dir, file_id_train, "_gateway_cellRef.rds")
    cell_state_df = readRDS(out_file)

    bk.dat = t(bk.dat)
    sc.dat = t(sc.dat)
    colnames(bk.dat) = colnames(sc.dat) 

    cell.type.labels = cell_state_df$cell_type
    cell.state.labels = cell.type.labels

    sc.stat <- plot.scRNA.outlier(
        input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
        cell.type.labels=cell.type.labels,
        species="hs", #currently only human(hs) and mouse(mm) annotations are supported
        return.raw=TRUE #return the data used for plotting.
    )

    bk.stat <- plot.bulk.outlier(
        bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID
        sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
        cell.type.labels=cell.type.labels,
        species="hs", #currently only human(hs) and mouse(mm) annotations are supported
        return.raw=TRUE
    )

    sc.dat.filtered <- cleanup.genes(input=sc.dat,
        input.type="count.matrix",
        species="hs",
        gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
        exp.cells=5)

    sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered,
        gene.type = "protein_coding")


    diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
        cell.type.labels=cell.type.labels,
        cell.state.labels=cell.state.labels,
        pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1cell.count.cutoff=50, # a numeric value to exclude cell state with numn.cores=1 #number of threads
    )

    sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
        stat=diff.exp.stat,
        pval.max=0.01,
        lfc.min=0.1)


    myPrism <- new.prism(
        reference=sc.dat.filtered.pc,
        mixture=bk.dat,
        input.type="count.matrix",
        cell.type.labels = cell.type.labels,
        cell.state.labels = cell.state.labels,
        key=NULL,
        outlier.cut=0.01,
        outlier.fraction=0.1,
    )

    bp.res <- run.prism(prism = myPrism, n.cores=1)

    out_file = paste0(out_dir, file_id_train, "_bp.res.rdata")    
    save(bp.res, file=out_file)

    Z.tumor <- get.exp(bp=bp.res, state.or.type="type")

    z_A549 = Z.tumor[1:dim(Z.tumor)[1], 1:dim(Z.tumor)[2],1]
    z_K562 = Z.tumor[1:dim(Z.tumor)[1], 1:dim(Z.tumor)[2],2]
    z_MCF7 = Z.tumor[1:dim(Z.tumor)[1], 1:dim(Z.tumor)[2],3]

    z_A549 = as.data.frame(z_A549)
    z_A549$cell_type = "A549"
    z_A549$samp_ID = 0:(dim(Z.tumor)[1]-1)

    z_K562 = as.data.frame(z_K562)
    z_K562$cell_type = "K562"
    z_K562$samp_ID = 0:(dim(Z.tumor)[1]-1)

    z_MCF7 = as.data.frame(z_MCF7)
    z_MCF7$cell_type = "MCF7"
    z_MCF7$samp_ID = 0:(dim(Z.tumor)[1]-1)

    dflist = list(z_A549, z_K562, z_MCF7)
    z_final = do.call("rbind",dflist)

    out_file = paste0(out_dir, file_id_train, "_bp_expr.tsv")    
    write.table(z_final, out_file, sep="\t", quote=F, row.names=F)

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


    res = format_bp_Z("A549", bp_res_Z)
    ct_z = res[[1]]
    ct_meta = res[[2]]


    res = format_bp_Z("K562", bp_res_Z)
    ct_z = rbind(ct_z, res[[1]])
    ct_meta = rbind(ct_meta, res[[2]])

    res = format_bp_Z("MCF7", bp_res_Z)
    ct_z = rbind(ct_z, res[[1]])
    ct_meta = rbind(ct_meta, res[[2]])


    outfile = paste0(out_dir, "bp_tocilizumab_Z_01.tsv")
    write.table(ct_z, outfile, sep="\t", quote=F, row.names=F)


    outfile = paste0(out_dir, "bp_tocilizumab_Z_meta_01.tsv")
    write.table(ct_meta, outfile, sep="\t", quote=F, row.names=F)

}


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


args = commandArgs(trailingOnly=TRUE)
file_id_train = args[1]
file_id_test = args[2]
in_dir = args[3]
out_dir = args[4]
num_samp = args[5]
ncores = as.numeric(args[6])

print("file_id_train:")
print(file_id_train)

print("file_id_test:")
print(file_id_test)

print("in_dir:")
print(in_dir)

print("out_dir:")
print(out_dir)

print("num_samp:")
print(num_samp)

print("ncores:")
print(ncores)

work_dir=${PWD}
in_dir=${work_dir}/../data/single_cell_data/cibersort_kang/
out_dir=${work_dir}/../results/single_cell_data/bp_kang/
file_id_train=bp
file_id_test=bp



format_bayes_prism(in_dir, out_dir, file_id_train, file_id_train)

run_bp(out_dir, file_id_train, ncores=1)





