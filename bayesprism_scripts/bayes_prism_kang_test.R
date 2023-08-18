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
    
    count_matr_t = t(count_matr)
    colnames(count_matr_t) = gene_ids
    rownames(count_matr_t) = cell_ids
    
    return(count_matr_t)
    
}

write_for_getway <- function(sc_ref, bulk_ref, out_dir, file_id_train){

    out_bulk_ref = t(bulk_ref)

    # make the cell reference matrix
    cell_id = c(1:nrow(sc_ref))
    cell_ref_df = data.frame(cell_id)
    cell_ref_df$cell_type = rownames(sc_ref)
    cell_ref_df$cell_subtype = rownames(sc_ref)
    cell_ref_df$tumor_state = rep(0, nrow(sc_ref))

    # make same cell id names
    rownames(sc_ref) = cell_id
    out_sc_ref = t(sc_ref)

    out_file = paste0(out_dir, file_id_train, "_gateway_bulk.tsv")
    write.table(out_bulk_ref, out_file, sep="\t", quote=F, row.names = T)
    out_file = paste0(out_dir, file_id_train, "_gateway_bulk.rds")
    saveRDS(as.data.frame(out_bulk_ref), out_file)

    out_file = paste0(out_dir, file_id_train, "_gateway_scRef.tsv")
    write.table(out_sc_ref, out_file, sep="\t", quote=F, row.names = T)
    out_file = paste0(out_dir, file_id_train, "_gateway_scRef.rds")
    saveRDS(as.data.frame(out_sc_ref), out_file)

    out_file = paste0(out_dir, file_id_train, "_gateway_cellRef.csv")
    write.table(cell_ref_df, out_file, sep=",", quote=F, row.names = T)
    out_file = paste0(out_dir, file_id_train, "_gateway_cellRef.rds")
    saveRDS(as.data.frame(cell_ref_df), out_file)


}

format_bayes_prism <- function(in_dir, out_dir, file_id_train, file_id_test, num_samp, ncores){
    sc_matr_file = paste0(in_dir, file_id_train, "_cybersort_sig.tsv.gz")
    sc_ref = read_and_format(sc_matr_file)
    
    bulk_matr_file = paste0(in_dir, file_id_test, "_cybersort_mix.tsv.gz")
    bulk_ref = read_and_format(bulk_matr_file)
    bulk_ref = data.frame(bulk_ref)
    rownames(bulk_ref) = paste0("samp_", rownames(bulk_ref))
    bulk_ref = as.matrix(bulk_ref)

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

    bp.res <- run.prism(prism = myPrism, n.cores=ncores)

    out_file = paste0(out_dir, "bp.res.rdata")    
    save(bp.res, file=out_file)

    Z.tumor <- get.exp(bp=bp.res, state.or.type="type")

    z_cd14 = Z.tumor[1:120, 1:11255,1]
    z_dc = Z.tumor[1:120, 1:11255,2]
    z_cd4_memT = Z.tumor[1:120, 1:11255,3]
    z_T_Act = Z.tumor[1:120, 1:11255,4]
    z_CD4_Naive_T = Z.tumor[1:120, 1:11255,5]
    z_CD8_T = Z.tumor[1:120, 1:11255,6]
    z_Mk = Z.tumor[1:120, 1:11255,7]
    z_B = Z.tumor[1:120, 1:11255,8]
    z_CD16_Mono = Z.tumor[1:120, 1:11255,9]
    z_NK = Z.tumor[1:120, 1:11255,10]

    z_cd14 = as.data.frame(z_cd14)
    z_cd14$cell_type = "CD14_Mono"
    z_cd14$samp_ID = 0:119

    z_dc = as.data.frame(z_dc)
    z_dc$cell_type = "DC"
    z_dc$samp_ID = 0:119

    z_cd4_memT = as.data.frame(z_cd4_memT)
    z_cd4_memT$cell_type = "CD4_Mem_T"
    z_cd4_memT$samp_ID = 0:119

    z_T_Act = as.data.frame(z_T_Act)
    z_T_Act$cell_type = "T_Act"
    z_T_Act$samp_ID = 0:119

    z_CD4_Naive_T = as.data.frame(z_CD4_Naive_T)
    z_CD4_Naive_T$cell_type = "CD4_Naive_T"
    z_CD4_Naive_T$samp_ID = 0:119

    z_CD8_T = as.data.frame(z_CD8_T)
    z_CD8_T$cell_type = "CD8_T"
    z_CD8_T$samp_ID = 0:119

    z_Mk = as.data.frame(z_Mk)
    z_Mk$cell_type = "Mk"
    z_Mk$samp_ID = 0:119

    z_B = as.data.frame(z_B)
    z_B$cell_type = "B"
    z_B$samp_ID = 0:119

    z_CD16_Mono = as.data.frame(z_CD16_Mono)
    z_CD16_Mono$cell_type = "CD16_Mono"
    z_CD16_Mono$samp_ID = 0:119

    z_NK = as.data.frame(z_NK)
    z_NK$cell_type = "NK"
    z_NK$samp_ID = 0:119

    dflist = list(z_cd14, z_dc, z_cd4_memT, z_T_Act, z_CD4_Naive_T, z_CD8_T, z_Mk, z_B, z_CD16_Mono, z_NK)
    z_final = do.call("rbind",dflist)

    out_file = paste0(out_dir, "bp_expr.tsv")    
    write.table(z_final, out_file, sep="\t", quote=F, row.names=F)

}


#!/usr/bin/env Rscript
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

format_bayes_prism(in_dir, out_dir, file_id_train, file_id_test, num_samp, ncores)
run_bp(bulk_file, scref_file, cellref_file, out_dir, ncores)

