#' SpaProxTarget
# 'Calculate the spatial proximity of each cell to the specified type (such as Tumor) and write back to the main table.
#' Based on the Euclidean distance, excluding itself, return:
#' -dist_to_nearest: The straight-line distance (µm) to the nearest target cell
#' -count_in_radius: Count of target cells within the radius (excluding itself)
#'
#' @param cell_data data.frame; It must contain four columns: x, y, cell_id and type
#' @param type Cell type column name (default "type")
#' @param target_type Target cell type, such as "Tumor" (default)
#' @param radius Circular search radius (µm, default 50)
# '
#' @returns Two new columns have been added to the original table: dist_to_nearest and count_in_radius
# '
#' @export
# '
#' @examples
#' target_cells <- SpaProxTarget(cell_data, type, "Tumor", 50)
SpaProxTarget <- function(cell_data,
                          type,
                          target_type = "Tumor",
                          radius = 50) {
  all_coords <- cell_data[, c("x", "y")]

  target_cells <- subset(cell_data, type == target_type)
  target_coords <- target_cells[, c("x", "y")]

  # If there are no target cells, return 0
  if (nrow(target_coords) == 0)
    return(rep(0, nrow(cell_data)))

  # Calculate the distance from each cell to the nearest target cell
  nn_dist <- nn2(target_coords, all_coords, k = 1)$nn.dists[, 1]

  # Calculate the number of target cells within the radius
  nn_count <- sapply(1:nrow(all_coords), function(i) {
    distances <- sqrt((target_coords$x - all_coords$x[i])^2 +
                        (target_coords$y - all_coords$y[i])^2)
    sum(distances <= radius)
  })
  cell_data$dist_to_nearest <- nn_dist
  cell_data$count_in_radius <- nn_count
  target_cells <- subset(cell_data, !(type == target_type))
  target_cells
}
