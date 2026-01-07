#' SpaComLeiden
#'
#' @param cell_data data.frame；It must contain four columns: x, y, cell_id and type
#' @param edges data.frame；It must contain three columns: from, to, and weight
#' @param resolution Resolution parameter (default 0.5
#' @param n_iter Number of iterations (default: 10
#'
#' @returns The igraph object has added the properties of community, degree, betweenness and closeness
#' @export
#' @examples
#' g <- SpaComLeiden(cell_data, edges, resolution = 0.5, n_iter = 10)
SpaComLeiden <- function(cell_data,
                         edges,
                         resolution = 0.5,
                         n_iter = 10) {
  g <- graph_from_data_frame(d = edges,
                             directed = FALSE,
                             vertices = cell_data)

  # Leiden zoning
  communities <- cluster_leiden(
    g,
    weights = E(g)$weight,
    objective_function = "modularity",
    resolution_parameter = resolution,
    n_iterations = n_iter
  )

  # Add topological attributes
  V(g)$community <- as.factor(membership(communities))
  V(g)$degree    <- degree(g)
  V(g)$betweenness <- betweenness(g)
  V(g)$closeness   <- closeness(g)

  g
}
