#' Plot Admixture-style plots of HTO proportions in droplets
#'
#' @param seurat_obj A Seurat object after HTODemux was run.
#' @param color_vec A named vector of custom colors, with HTOs as names.
#'
#' @return A ggplot.
#' @export
#'
#' @examples
#' mycolors <- rainbow(4)
#' names(mycolors) <- c("black", "green", "blue", "tomato3")
#' plot_admix(pbmc.3k, colors_vec = rainbow(8))
plot_admix <- function(seurat_obj, colors_vec) {
    
    meta_df <- seurat_obj@meta.data |>
        as_tibble(rownames = "barcode") |>
        select(barcode, stim = HTO_maxID, hto_class = HTO_classification.global)
  
    hto_df <- seurat_obj@assays$HTO@counts |>
        as_tibble(rownames = "hto") |>
        pivot_longer(-hto, names_to = "barcode") |>
        left_join(meta_df)
    
    hto_top <- hto_df |>
        group_by(barcode, hto_class) |>
        slice_max(n = 1, order_by = value) |>
        select(barcode, hto_class, top_hto = hto) |>
        ungroup()
    
    hto_plot_df <- hto_df |>
        left_join(hto_top, by = c("barcode", "hto_class")) |>
        mutate_at(vars(top_hto, stim), ~factor(., levels = names(colors_vec)))
    
    ggplot(hto_plot_df, aes(reorder_within(barcode, by = value, within = stim), value)) +
    geom_col(aes(fill = hto), position = "fill", width = 1.01, show.legend = FALSE) +
    scale_fill_manual(values = colors_vec) +
    facet_grid(~stim, scales = "free_x", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.text = element_blank()) +
    labs(x = NULL, y = NULL)
}
