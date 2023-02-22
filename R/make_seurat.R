#' Make Seurat object 
#'
#' @param cellranger_path Path to the CellRanger filtered feature matrix directory.
#' @param project_id A name for your project.
#' @param hto_names A vector for naming each hashtag, in the same order as in the feature matrix.
#' @param mito_ids A vector of ENSEMBL gene IDs for mitochondrial genes.
#' @param ribo_ids A vector of ENSEMBL gene IDs for ribosomal protein genes.
#'
#' @return A Seurat Object
#'
#' @export
#'
#' @examples
#' my_obj <- make_seurat("filtered_feature_bc_matrix", "PBMCs", c("condition1", "condition2"), mito_genes_ids, ribo_genes_ids)

make_seurat <- function(cellranger_path, project_id, hto_names, mito_ids, ribo_ids) {

    data10x <- Read10X(cellranger_path, gene.column = 1)

    antibody_mtx <- data10x[["Antibody Capture"]] |>
        {function(x) x[!grepl("^Hashtag", rownames(x)), ]}()

    rownames(antibody_mtx) <- sub("_prot$", "", rownames(antibody_mtx))

    hashtags_mtx <- data10x[["Antibody Capture"]] |>
        {function(x) x[grepl("^Hashtag", rownames(x)), ]}()

    rownames(hashtags_mtx) <- setNames(hto_names[rownames(hashtags_mtx)], NULL)

    bcells <- CreateSeuratObject(counts = data10x[["Gene Expression"]],
                                 project = project_id)

    bcells[["ADT"]] <- CreateAssayObject(counts = antibody_mtx)
    bcells[["HTO"]] <- CreateAssayObject(counts = hashtags_mtx)

    bcells <- bcells |>
        NormalizeData(normalization.method = "LogNormalize", margin = 1) |>
        NormalizeData(assay = "HTO", normalization.method = "CLR", margin = 2) |>
        NormalizeData(assay = "ADT", normalization.method = "CLR", margin = 2)

    bcells[["percent_mt"]] <- PercentageFeatureSet(bcells, features = mito_ids)
    bcells[["percent_ribo"]] <- PercentageFeatureSet(bcells, features = ribo_ids)

    bcells
}
