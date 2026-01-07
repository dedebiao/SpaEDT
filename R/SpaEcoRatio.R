#' SpaEcoRatio
# 'Calculate the "K-nearest neighbor ecological ratio" for spatial single-cell data.
# 'Take k nearest neighbors (excluding the cell itself) for each cell and calculate the proportion of each cell type.
# '
#' Core Logic:
#' 1. KNN finds neighbors (excluding itself);
#' 2. Divide the type count by k to get the proportion;
#' 3. Write back to the main table without destroying the original columns.
# '
#' @param cell_data data.frame; It must contain x, y (unit: µm) and any other column
#' @param type Cell type column name (default "type")
#' @param x, y coordinate column names (default "x", "y")
#' @param k Number of nearest neighbors (default 50, excluding itself)
# '
#' @returns The original table adds a new column "{type}_proportion" (0-1 consecutive)
# '
#' @export
#'
#' @examples
#' cell_data <- SpaEcoRatio(cell_data,x = "x",y = "y",type = "type",k=50)
SpaEcoRatio <- function(cell_data,
                        x = "x",
                        y = "y",
                        type = "type",
                        k = 50) {
  # 1. Find the nearest k-neighbor (excluding yourself) in a k-d tree
  knn_out <- RANN::nn2(data = cell_data[c("x", "y")], k = k + 1)
  nb_idx  <- knn_out$nn.idx[, -1]

  # 2. Count the proportion of neighbor types of each cell
  prop_lst <- lapply(seq_len(nrow(cell_data)), function(i) {
    nei_anno <- cell_data[[type]][nb_idx[i, ]]
    prop_tbl <- table(factor(nei_anno, levels = unique(cell_data[[type]]))) / k
    as.list(prop_tbl)
  })

  # 3. Bind back to the main table
  prop_cell_data <- bind_rows(prop_lst) %>% rename_with( ~ paste0(., "_proportion"))
  cell_data  <- bind_cols(cell_data, prop_cell_data)
}
