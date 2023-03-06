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

make_seurat <- function(cellranger_path, project_id, hto_names = NULL, mito_ids, ribo_ids) {

    data10x <- Read10X(cellranger_path, gene.column = 1)

    antibody_mtx <- data10x[["Antibody Capture"]] |>
        {function(x) x[!grepl("^Hashtag", rownames(x)), ]}()


    if (!is.null(hto_names)) { 

	hashtags_mtx <- data10x[["Antibody Capture"]] |>
	    {function(x) x[grepl("^Hashtag", rownames(x)), ]}()

	rownames(hashtags_mtx) <- setNames(hto_names[rownames(hashtags_mtx)], NULL)
    }

    seuratobj <- 
	CreateSeuratObject(counts = data10x[["Gene Expression"]], project = project_id) |>
        NormalizeData(assay = "RNA", normalization.method = "LogNormalize")

    if ( nrow(antibody_mtx) > 0 ) {
	
	rownames(antibody_mtx) <- sub("_prot$", "", rownames(antibody_mtx))
	
	seuratobj[["ADT"]] <- CreateAssayObject(counts = antibody_mtx)
    
	seuratobj <- NormalizeData(seuratobj, assay = "ADT", normalization.method = "CLR", margin = 2)
    }
    
    if (!is.null(hto_names)) {
	
	seuratobj[["HTO"]] <- CreateAssayObject(counts = hashtags_mtx)

	seuratobj <- NormalizeData(seuratobj, assay = "HTO", normalization.method = "CLR", margin = 2)
    }

    seuratobj[["percent_mt"]] <- PercentageFeatureSet(seuratobj, features = mito_ids)
    seuratobj[["percent_ribo"]] <- PercentageFeatureSet(seuratobj, features = ribo_ids)

    seuratobj
}
