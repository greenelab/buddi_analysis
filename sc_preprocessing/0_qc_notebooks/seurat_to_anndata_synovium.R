library(Seurat)
library(SeuratDisk)

curr_dir = getwd()


col_interest = c("cell", "sample", "cluster_number", "cluster_name")

# the metadata files are stored as independent seurat objects
# we will get the metadata, used it to subset the full raw data matrix
# then combine into a single seurat object

get_raw_counts <- function(metadata_file, raw_data_file, col_interest, cell_type){

    full_raw_data = readRDS(raw_data_file)

    curr_cell = readRDS(metadata_file)

    # subset data to only current cell type
    curr_cell_raw = full_raw_data[,curr_cell$meta_data$cell]

    # get metadata and format
    meta_data_df = curr_cell$meta_data
    rownames(meta_data_df) = curr_cell$meta_data$cell

    # remove columns that are not the same across all cell types
    meta_data_df = meta_data_df[,col_interest]
    meta_data_df$cell_type = cell_type

    # make SeuratObject
    curr_seurat = CreateSeuratObject(
        curr_cell_raw,
        assay = "RNA",
        meta.data = meta_data_df,
        project = "SeuratProject",
        min.cells = 0,
        min.features = 0
    )

    # return
    return(curr_seurat)
}

raw_data_file = paste0(curr_dir, "/data/single_cell_data/synovium/raw_mRNA_count_matrix.rds")

nk_metadata_file = paste0(curr_dir, "/data/single_cell_data/synovium/NK_reference.rds")
nk_seurat = get_raw_counts(nk_metadata_file, raw_data_file, col_interest, cell_type = "NK")


metadata_file = paste0(curr_dir, "/data/single_cell_data/synovium/Bcell_reference.rds")
b_seurat = get_raw_counts(metadata_file, raw_data_file, col_interest, cell_type = "B")

full_seurat = merge(nk_seurat, b_seurat)
rm(list=c("nk_seurat", "b_seurat" ))


metadata_file = paste0(curr_dir, "/data/single_cell_data/synovium/endothelial_reference.rds")
endothelial_seurat = get_raw_counts(metadata_file, raw_data_file, col_interest, cell_type = "Endothelial")

full_seurat = merge(full_seurat, endothelial_seurat)
rm(list=c("endothelial_seurat" ))


metadata_file = paste0(curr_dir, "/data/single_cell_data/synovium/fibroblast_reference.rds")
fib_seurat = get_raw_counts(metadata_file, raw_data_file, col_interest, cell_type = "Fibroblast")

full_seurat = merge(full_seurat, fib_seurat)
rm(list=c("fib_seurat" ))


metadata_file = paste0(curr_dir, "/data/single_cell_data/synovium/myeloid_reference.rds")
myeloid_seurat = get_raw_counts(metadata_file, raw_data_file, col_interest, cell_type = "Myeloid")

full_seurat = merge(full_seurat, myeloid_seurat)
rm(list=c("myeloid_seurat" ))


metadata_file = paste0(curr_dir, "/data/single_cell_data/synovium/t_reference.rds")
t_seurat = get_raw_counts(metadata_file, raw_data_file, col_interest, cell_type = "T")

full_seurat = merge(full_seurat, t_seurat)
rm(list=c("t_seurat" ))

# convert and write to disc
seurat_filename = paste0(curr_dir, "/data/single_cell_data/synovium/synovium_raw_reference.h5Seurat")
SaveH5Seurat(full_seurat, filename = seurat_filename)
Convert(seurat_filename, dest = "h5ad")


