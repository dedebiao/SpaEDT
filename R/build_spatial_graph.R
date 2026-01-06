#' build_spatial_graph
#'
#' @param cell_data
#' @param cell_id
#' @param k
#' @param max_dist
#'
#' @returns
#' @export
#'
#' @examples
build_spatial_graph <- function(cell_data,
                                cell_id = 'cell_id',
                                k = 6,
                                max_dist = 100) {
  coords <- cell_data %>% dplyr::select(x, y) %>% as.matrix()

  knn_result <- get.knn(coords, k = k)

  edges <- data.frame()
  for (i in 1:nrow(coords)) {
    neighbors <- knn_result$nn.index[i, ]
    distances <- knn_result$nn.dist[i, ]

    valid <- which(distances <= max_dist)
    if (length(valid) > 0) {
      new_edges <- data.frame(
        from = cell_data$cell_id[i],
        to = cell_data$cell_id[neighbors[valid]],
        weight = distances[valid],
        type = paste0(cell_data$type[i], "-", cell_data$type[neighbors[valid]])
      )
      edges <- rbind(edges, new_edges)
    }
  }

  g <- graph_from_data_frame(d = edges,
                             directed = FALSE,
                             vertices = cell_data)

  V(g)$degree <- degree(g)
  V(g)$betweenness <- betweenness(g)
  V(g)$closeness <- closeness(g)

  list(graph = g, edges = edges)
}
