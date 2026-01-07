#' SpaBarrierQ
# 'Quantify the spatial Barrier Effect of stromal cells on tumor-immune contact.
#For each immune cell, calculate:
#' -Blocking count: The number of stromal cells (excluding themselves) that are simultaneously located within the "radius around immune cells" and the "radius around nearest tumor cells".
# '
#' Principle:
#' - Take the nearest tumor cell as the anchor point and expand radius_factor × tumor distance outward as the search radius;
#Count the number of stromal cells within the intersection of these two radii;
#The higher the value, the stronger the matrix barrier, suggesting immune rejection or physical obstruction.
# '
#' @param cell_data data.frame; It must contain four columns: x, y, cell_id and type
#' @param tumor_type Anchor cell type, default "Tumor"
#' @param immune_type Blocked cell type, default "Immune"
#' @param stroma_type Barrier cell type, default "Stroma"
#' @param radius_factor search radius multiple (default 1.5), ≥1
# '
#' @returns The original table adds the column stroma_blocker (integer ≥ 0)
# '
#' @export
# '
#' @examples
#' cell_data <- SpaBarrierQ(cell_data,type,tumor_type = "Tumor",mmune_type = "Immune",stroma_type = "Stroma")
SpaBarrierQ <- function(cell_data,
                        type,
                        tumor_type = "Tumor",
                        immune_type = "Immune",
                        stroma_type = "Stroma",
                        radius_factor = 1.5) {
  tumor_cells <- cell_data %>% filter(type == tumor_type)
  immune_cells <- cell_data %>% filter(type == immune_type)
  stroma_cells <- cell_data %>% filter(type == stroma_type)

  tumor_tree <- nn2(tumor_cells[, c("x", "y")], k = 1)
  stroma_tree <- nn2(stroma_cells[, c("x", "y")], k = nrow(stroma_cells))

  blockers <- rep(NA, nrow(immune_cells))

  if (nrow(immune_cells) > 0 &&
      nrow(tumor_cells) > 0 && nrow(stroma_cells) > 0) {
    for (i in 1:nrow(immune_cells)) {
      immune_loc <- immune_cells[i, c("x", "y")]

      # Find the nearest tumor cell
      tumor_dist <- min(sqrt((tumor_cells$x - immune_loc$x)^2 +
                               (tumor_cells$y - immune_loc$y)^2
      ))

      # Calculate the search radius
      search_radius <- tumor_dist * radius_factor

      # Search for stromal cells around immune cells
      stroma_near_immune <- stroma_cells %>%
        mutate(dist = sqrt((x - immune_loc$x)^2 + (y - immune_loc$y)^2)) %>%
        filter(dist <= search_radius) %>%
        pull(cell_id)

      # Search for stromal cells around tumor cells
      tumor_near_immune <- tumor_cells %>%
        mutate(dist = sqrt((x - immune_loc$x)^2 + (y - immune_loc$y)^2)) %>%
        filter(dist <= tumor_dist) %>%
        slice(which.min(dist))

      if (nrow(tumor_near_immune) > 0) {
        stroma_near_tumor <- stroma_cells %>%
          mutate(dist = sqrt((x - tumor_near_immune$x)^2 +
                               (y - tumor_near_immune$y)^2
          )) %>%
          filter(dist <= search_radius) %>%
          pull(cell_id)

        # Calculate the blocked stromal cells
        blocker_cells <- intersect(stroma_near_immune, stroma_near_tumor)
        blockers[i] <- length(blocker_cells)
      } else {
        blockers[i] <- 0
      }
    }
  }

  result <- data.frame(cell_id = immune_cells$cell_id, stroma_blocker = blockers)

  cell_data <- left_join(cell_data, result, by = "cell_id")
  cell_data$stroma_blocker[is.na(cell_data$stroma_blocker)] <- 0

  cell_data
}
