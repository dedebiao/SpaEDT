#' SpaGeoNet
# '
#' Construct Spatial Geometry Network (Spatial Geometry Network).
#' Generate an undirected, weighted adjacency edge list based on the KNN + distance threshold
# '
#' Core Logic:
#'1. Take k nearest neighbors (excluding the cell itself) for each cell;
# 2. Only retain edges that are ≤ max_dist µm;
# '
#' @param cell_data data.frame; It must contain four columns: x, y, cell_id, and type (unit: µm).
#' @param cell_id Cell Unique ID column name (default "cell_id")
#' @param k Number of nearest neighbors (default 50, excluding itself)
#' @param max_dist Maximum side length (µm, default 100)
# '
#' @returns data.frame: from, to, weight, type four columns, directly graph_from_data_frame()
# '
#' @export
# '
#' @examples
#' edges <- SpaGeoNet(cell_data, cell_id='cell_id', k = 6, max_dist = 50)
SpaGeoNet <- function(cell_data,
                      cell_id = 'cell_id',
                      k = 50,
                      max_dist = 100) {
  coords <- cell_data %>% dplyr::select(x, y) %>% as.matrix()

  knn_result <- FNN::get.knn(coords, k = k)

  # Create an adjacency list
  edges <- data.frame()
  for (i in 1:nrow(coords)) {
    neighbors <- knn_result$nn.index[i, ]
    distances <- knn_result$nn.dist[i, ]

    # Apply the distance threshold
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

  edges
}
