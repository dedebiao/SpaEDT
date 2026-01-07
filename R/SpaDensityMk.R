#' SpaDensityMk
# 'Quickly generate physical density maps for spatial single-cell data.
# 'Using a circular window as the search unit, return the neighbor count and physical density (cells/mm²) excluding itself.
# '
#' Core Logic:
#' 1. Count of neighbors within the radius that do not include itself;
#' 2. Divide by the area of the circle (mm²) to obtain the physical density;
# '
#' @param cell_data data.frame; It must contain x, y (unit: µm) and any other column
#' @param radius Circular window Radius (µm, default 50)
# '
#' @returns The original table adds a new column density_mm2 (cells/mm²)
# '
#' @export
# '
#' @examples
#' cell_data <- SpaDensityMk(cell_data, radius = 20)
SpaDensityMk <- function(cell_data, radius = 50) {
  radius_mm <- radius / 1000                # mm
  area <- pi * radius_mm^2                  # mm²
  coords <- cell_data[, c("x", "y")]
  nn <- RANN::nn2(coords, coords, searchtype = "radius", radius = radius)
  nei_cnt <- sapply(1:nrow(coords), function(i)
    sum(nn$nn.idx[i, ] > 0) - 1)
  density_mm2 <- nei_cnt / area
  cell_data$density_mm2 <- density_mm2
  cell_data
}
