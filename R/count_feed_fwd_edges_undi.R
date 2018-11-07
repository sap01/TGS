#' Count the number of feed-forward edges in a given undirected network.
#'
#' @param undi.net.adj.matrix adjacency matrix of an undirected network
#'
#' @return count of the number of feed-forward edges
#' @examples
#' CountFeedFwdEdgesUndi(matrix(c(0,1,0,1,0,1,0,0,0),nrow=3))
#'
#' @export
CountFeedFwdEdgesUndi <- function(undi.net.adj.matrix) {
  if(!base::is.matrix(undi.net.adj.matrix))
  {
    base::stop("Error in CountFeedFwdEdgesUndi undi.net.adj.matrix is not a matrix")
  }
  num.nodes <- base::nrow(undi.net.adj.matrix)

  count <- 0

  for (node.idx.1 in 3:num.nodes) {
    for (node.idx.2 in 2:(node.idx.1 - 1)) {
      for (node.idx.3 in 1:(node.idx.2 - 1)) {

        if ((undi.net.adj.matrix[node.idx.1, node.idx.2] == 1) &
          (undi.net.adj.matrix[node.idx.2, node.idx.3] == 1) &
          (undi.net.adj.matrix[node.idx.1, node.idx.3] == 1)) {

          count <- (count + 1)
        }
      }
      base::rm(node.idx.3)
    }
    base::rm(node.idx.2)
  }
  base::rm(node.idx.1)

  return(count)

}
