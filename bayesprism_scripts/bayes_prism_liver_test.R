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
        species="mm",
        gene.group=c( "Rb","Mrp","other_Rb","chrM","chrX","chrY") ,
        exp.cells=5)

    # skipping this step because it doesn't look like its working for mouse
    #sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered,
    #    gene.type = "protein_coding")


    diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
        cell.type.labels=cell.type.labels,
        cell.state.labels=cell.state.labels,
        pseudo.count=0.1, 
    )

    sc.dat.filtered.sig <- select.marker (sc.dat=sc.dat.filtered,
        stat=diff.exp.stat,
        pval.max=0.01,
        lfc.min=0.1)


    myPrism <- new.prism(
        reference=sc.dat.filtered,
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

    hepatocyte = Z.tumor[1:49, 1:11052,1]
    hepatic_sinusoid = Z.tumor[1:49, 1:11052,2]
    kupffer = Z.tumor[1:49, 1:11052,3]
    hepatic_stellate = Z.tumor[1:49, 1:11052,4]
    NK = Z.tumor[1:49, 1:11052,5]
    plasmacytoid_dendritic_cell = Z.tumor[1:49, 1:11052,6]
    b_cell = Z.tumor[1:49, 1:11052,7]
    myeloid_leukocyte = Z.tumor[1:49, 1:11052,8]

    hepatocyte = as.data.frame(hepatocyte)
    hepatocyte$cell_type = "hepatocyte"
    hepatocyte$sample_id = row.names(hepatocyte)

    hepatic_sinusoid = as.data.frame(hepatic_sinusoid)
    hepatic_sinusoid$cell_type = "hepatic_sinusoid"
    hepatic_sinusoid$sample_id = row.names(hepatic_sinusoid)


    kupffer = as.data.frame(kupffer)
    kupffer$cell_type = "kupffer"
    kupffer$sample_id = row.names(kupffer)


    hepatic_stellate = as.data.frame(hepatic_stellate)
    hepatic_stellate$cell_type = "hepatic_stellate"
    hepatic_stellate$sample_id = row.names(hepatic_stellate)


    NK = as.data.frame(NK)
    NK$cell_type = "NK"
    NK$sample_id = row.names(NK)


    plasmacytoid_dendritic_cell = as.data.frame(plasmacytoid_dendritic_cell)
    plasmacytoid_dendritic_cell$cell_type = "plasmacytoid_dendritic_cell"
    plasmacytoid_dendritic_cell$sample_id = row.names(plasmacytoid_dendritic_cell)


    b_cell = as.data.frame(b_cell)
    b_cell$cell_type = "b_cell"
    b_cell$sample_id = row.names(b_cell)


    myeloid_leukocyte = as.data.frame(myeloid_leukocyte)
    myeloid_leukocyte$cell_type = "myeloid_leukocyte"
    myeloid_leukocyte$sample_id = row.names(myeloid_leukocyte)



    dflist = list(hepatocyte, hepatic_sinusoid, kupffer, hepatic_stellate, NK, plasmacytoid_dendritic_cell, b_cell, myeloid_leukocyte)
    z_final = do.call("rbind",dflist)

    z_final$sample_id = sub("samp_", "", z_final$sample_id)

    out_file = paste0(out_dir, "bp_expr.tsv")    
    write.table(z_final, out_file, sep="\t", quote=F, row.names=F)

}


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

in_dir = args[1]
out_dir = args[2]
file_id_train = args[3]
file_id_test = args[4]
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


run_bp(out_dir, file_id_train, ncores)

