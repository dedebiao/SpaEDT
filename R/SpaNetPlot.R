#' SpaNetPlot Draw the spatial community network
#'
#' @param g igraph object; It must contain attributes such as x, y, name, community, and degree
#' @param node_alpha node transparency (default 0.7)
#' @param title Graph Title (default: "Spatial Community Cluster")
#'
#' @return ggplot object
#' @export
#' @examples
#' SpaNetPlot(g, node_alpha = 0.7)
SpaNetPlot <- function(g,
                       node_alpha = 0.7,
                       title = "Spatial Community Clsuter") {
  # 2. 节点数据框
  node_df <- data.frame(
    cell_id = V(g)$name,
    x       = V(g)$x,
    y       = V(g)$y,
    degree  = degree(g),
    community = V(g)$community
  )

  ggplot() +
    geom_point(data = node_df,
               aes(
                 x = x,
                 y = y,
                 color = community,
                 size = degree
               ),
               alpha = node_alpha) +
    scale_size_continuous(range = c(1, 5)) +
    theme_void() +
    labs(title = title, x = "", y = "") +
    guides(size = guide_legend(title = "Degree"),
           color = guide_legend(title = "Community"))
}
