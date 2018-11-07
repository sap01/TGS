#' Returns all the nodes reachable from the given node in the directed adjacency matrix
#'
#' Given a directed network adjacency matrix and a node name, returns names of all the nodes
#' reachable from the given node. In the given directed network adjacency matrix,
#' rows = source nodes, cols = tgt nodes. (i, j)th cell = 1 implies a directed edge from i
#' to j; = 0 implies no directed edge from i to j. Row names of the given matrix must be
#' the names of the nodes.
#'
#' @param di.net.adj.matrix The directed network adjacency matrix
#' @param src.node.name The node whose reachable nodes need to be returned
#'
#' @return A vector containing all the nodes reachable from the given node
#' @examples
#' x = matrix(c(1,0,0,0,1,1,0,1,1),nrow=3)
#' rownames(x) <- c('A','B','C')
#' colnames(x) <- c('A','B','C')
#' reachable.nodes(x,'A')
#' @export
reachable.nodes <- function(di.net.adj.matrix, src.node.name)
{
  if(!base::is.matrix(di.net.adj.matrix))
  {
    base::stop("Error in reachable.nodes. di.net.adj.matrix is not a matrix")
  }
  reachable.node.names <- base::c()

  ## If the given node has at least one target node
  if (base::length(base::which(di.net.adj.matrix[src.node.name, ] == 1)) > 0)
  {
    reachable.node.names <- base::names(which(di.net.adj.matrix[src.node.name, ] == 1))

    ## Avoid self loop
    to.traverse.node.names <- base::setdiff(reachable.node.names, src.node.name)

    while (base::length(to.traverse.node.names) != 0)
    {
      curr.src.node.name <- to.traverse.node.names[1]

      if (base::length(base::which(di.net.adj.matrix[curr.src.node.name, ] == 1)) > 0)
      {
        new.reachable.node.names <- base::names(which(di.net.adj.matrix[curr.src.node.name, ] == 1))

        ## Avoid directed cycle
        new.reachable.node.names <- base::setdiff(new.reachable.node.names, reachable.node.names)

        reachable.node.names <- base::union(reachable.node.names, new.reachable.node.names)

        to.traverse.node.names <- base::union(to.traverse.node.names, new.reachable.node.names)
      }

      to.traverse.node.names <- base::setdiff(to.traverse.node.names, curr.src.node.name)
    }
  }

  return(reachable.node.names)
}
