library(Seurat)
library(SeuratDisk)
library(monocle3)
library(sceasy)
library(reticulate)

# this is using python  3.8.11
# you have to use your base installation of python on OS X
# this is because when python is built it needs the --enable-sharing
# flag to be set, but this is not set in conda  or virtualenvs on OS X
# this means that in your base installation of python, you need AnnData

use_python("/opt/homebrew/Caskroom/miniconda/base/bin/python/")

curr_dir = getwd()



# we have a three monocle 3 objects that we will write out the metadata tables for
# and write out the expression file as an Seurat --> AnnData object
# so we can read it into Python easily

# the files are divided by cell line

get_raw_counts <- function(raw_data_file){

    monocle_df = readRDS(raw_data_file)

    metadata_df = colData(df)

    gene_df = rowData(df)

    # make SeuratObject
    seurat_df = CreateSeuratObject(exprs(monocle_df), project = "SeuratProject", assay = "RNA",
                                    min.cells = 0, min.features = 0,
                                    names.delim = "-")


    return(list(seurat_df, metadata_df, gene_df))
}


raw_data_file = paste0(curr_dir, "/data/single_cell_data/sciplex/GSM4150378_sciPlex3_A549_24hrs.rds")
monocle_df = readRDS(raw_data_file)

outname = paste0(curr_dir, "/data/single_cell_data/sciplex/GSM4150378_sciPlex3_A549_24hrs.h5ad")
sceasy::convertFormat(monocle_df, from="sce", to="anndata",
                       outFile=outname)




raw_data_file = paste0(curr_dir, "/data/single_cell_data/sciplex/GSM4150378_sciPlex3_K562_24hrs.rds")
monocle_df = readRDS(raw_data_file)

outname = paste0(curr_dir, "/data/single_cell_data/sciplex/GSM4150378_sciPlex3_K562_24hrs.h5ad")
sceasy::convertFormat(monocle_df, from="sce", to="anndata",
                       outFile=outname)



raw_data_file = paste0(curr_dir, "/data/single_cell_data/sciplex/GSM4150378_sciPlex3_MCF7_24hrs.rds")
monocle_df = readRDS(raw_data_file)

outname = paste0(curr_dir, "/data/single_cell_data/sciplex/GSM4150378_sciPlex3_MCF7_24hrs.h5ad")
sceasy::convertFormat(monocle_df, from="sce", to="anndata",
                       outFile=outname)




res = get_raw_counts(raw_data_file)
a549_expr = res[[1]]
a549_col = res[[2]]
a549_gene = res[[3]]

raw_data_file = paste0(curr_dir, "/data/single_cell_data/sciplex/GSM4150378_sciPlex3_K562_24hrs.rds")
res = get_raw_counts(raw_data_file)
k562_expr = res[[1]]
k562_col = res[[2]]
k562_gene = res[[3]]

full_seurat = merge(a549_expr, k562_expr)
rm(list=c("a549_expr", "k562_expr" ))


raw_data_file = paste0(curr_dir, "/data/single_cell_data/sciplex/GSM4150378_sciPlex3_MCF7_24hrs.rds")
res = get_raw_counts(raw_data_file)
mcf7_expr = res[[1]]
mcf7_col = res[[2]]
mcf7_gene = res[[3]]

full_seurat = merge(full_seurat, mcf7_expr)
rm(list=c("mcf7_expr"))


# convert and write to disc
seurat_filename = paste0(curr_dir, "/data/single_cell_data/sciplex/sciplex_raw_reference.h5Seurat")
SaveH5Seurat(monocle_df, filename = seurat_filename)
Convert(seurat_filename, dest = "h5ad")


# make sure gene tables are the same
sum(a549_gene[,1] != k562_gene[,1]) == 0
sum(mcf7_gene[,1] != k562_gene[,1]) == 0

# make sure metadata columns are the same
sum(colnames(a549_col) != colnames(k562_col)) == 0
sum(colnames(mcf7_col) != colnames(k562_col)) == 0

metadata_total = rbind(a549_col, k562_col)
metadata_total = rbind(metadata_total, mcf7_col)

# write out the metadata
cellmeta_filename = paste0(curr_dir, "/data/single_cell_data/sciplex/cell_meta.tsv")
write.table(metadata_total, cellmeta_filename, sep="\t", quote=F, row.names=F)

genemeta_filename = paste0(curr_dir, "/data/single_cell_data/sciplex/gene_meta.tsv")
write.table(a549_gene, genemeta_filename, sep="\t", quote=F, row.names=F)

