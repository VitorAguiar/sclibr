#' Plot marker genes
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_df A data.frame such as returned by Seurat::FindAllMarkers with added gene_id and gene_name columns.
#'
#' @return A dotplot of the top 10 marker genes per identity.
#' 
#' @export
#'
#' @examples
#'
plot_markers <- function(seurat_obj, cluster_df) {
  
    top_markers <- cluster_df |>
        as_tibble() |>
        group_by(cluster) |>
        top_n(10, avg_log2FC) |>
        ungroup() |>
        select(cluster, gene_name, gene_id = gene, avg_log2FC)
  
    cluster_cell <- Idents(seurat_obj) |>
        enframe(name = "barcode", value = "cluster_cell") |>
        mutate(cluster_cell = factor(cluster_cell, levels = levels(top_markers$cluster)))
  
    cell_gene_expr <- seurat_obj@assays$RNA@data |>
        {function(x) x[unique(top_markers$gene_id), ]}() |>
        as_tibble(rownames = "gene_id") |>
        pivot_longer(-gene_id, names_to = "barcode", values_to = "logexpr")
  
    top_marker_expr <- cell_gene_expr |>
        left_join(cluster_cell, by = "barcode") |>
        inner_join(top_markers, by = "gene_id", multiple = "all")
  
    top_marker_perc_exp <- top_marker_expr |>
        group_by(cluster = cluster_cell, gene_name) |>
        summarise(prop_expr = mean(logexpr > 0)) |>
        ungroup()
  
    top_marker_avg_exp <- top_marker_expr |>
        filter(logexpr > 0) |>
        group_by(cluster = cluster_cell, gene_name) |>
        summarise(scaled_expr = mean(logexpr)) |>
        ungroup()
  
    top_marker_summary <- left_join(top_marker_perc_exp, top_marker_avg_exp) |>
        left_join(select(top_markers, cluster_top = cluster, gene_name, avg_log2FC), 
                  multiple = "all") |>
        mutate_at(vars(cluster, cluster_top), factor) |>
        arrange(cluster_top, avg_log2FC) |>
        mutate(gene_name = fct_inorder(gene_name)) 
  
    ggplot(top_marker_summary, aes(cluster, gene_name)) +
        geom_point(aes(size = prop_expr, fill = scaled_expr), 
                   color = "black", shape = 21) +
        scale_size(range = c(0.1, 4), labels = scales::percent) +
        scale_fill_viridis_c(option = "magma") +
        facet_wrap(~cluster_top, scales = "free", ncol = 3) +
        theme_bw() +
        theme(axis.line = element_blank(),
              axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 8),
              panel.grid = element_line(color = "grey96"),
              panel.border = element_blank(),
              legend.position = "top") +
        labs(x = NULL, y = NULL, 
             fill = "Scaled\nExpression", 
             size = "% of cells") +
        guides(fill = guide_colorbar(barheight = .5))
}

