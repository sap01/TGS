#' Generates True net adjacency matrix and save as an R object
#'
#' @param input.file full path of the input file containing adjacency list
#' @param output.file full path where output adjacency matrix is to be stored
#' @param num.nodes the number of nodes
#'
#' @import utils
#'
#' @return also returns the adjacency matrix
#'
#' @export
GenTrueAdjMatrix <- function(input.file, output.file, num.nodes)
{
  node.names <- base::c()
  for (node.idx in 1:num.nodes)
  {
    new.node.name <- base::paste('G', as.character(node.idx), sep = '')
    node.names <- base::c(node.names, new.node.name)
  }
  ## End: Specify input params

  # Param 'colClasses' must be given. Otherwise, the node names get converted into factors.
  true.net.adj.list <- utils::read.table(input.file, header = FALSE, sep = '\t',
                                  col.names = base::c('src.node', 'tgt.node', 'isEdge'),
                                  colClasses = base::c('character', 'character', 'integer'))

  ## Begin: Convert adjacency list to adjacency matrix
  true.net.adj.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                                dimnames = base::list(node.names, node.names))
  for (row.idx in 1:nrow(true.net.adj.list))
  {
    if (true.net.adj.list[row.idx, 'isEdge'] == 1)
    {
      src.node <- true.net.adj.list[row.idx, 'src.node']
      tgt.node <- true.net.adj.list[row.idx, 'tgt.node']
      true.net.adj.matrix[src.node, tgt.node] <- 1
    }
  }

  ## Check number of edges
  # length(true.net.adj.matrix[true.net.adj.matrix == 1])

  ## End: Convert adjacency list to adjacency matrix

  base::save(true.net.adj.matrix, file = output.file)

  return(true.net.adj.matrix)
}
