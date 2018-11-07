#' Checks if 'di.net.adj.matrix' = 'cmi.net.adj.matrix'
#'
#' @param di.net.adj.matrix First adjacency metric
#' @param cmi.net.adj.matrix Second adjacency metric
#' @param num.node Number of nodes in the metrices
#'
#' @return Returns 1 if 'di.net.adj.matrix' = 'cmi.net.adj.matrix' else returns 0
#' @examples
#' CompareNet(matrix(c(0,0,0,0),nrow=2,ncol=2),
#' + matrix(c(0,0,0,0),nrow=2,ncol=2),
#' + 2)
#'
#' CompareNet(matrix(c(0,0,1,0),nrow=2,ncol=2),
#' + matrix(c(0,0,0,0),nrow=2,ncol=2),
#' + 2)
#' @export
CompareNet <-function(di.net.adj.matrix,cmi.net.adj.matrix,num.node)
{
  if(!base::is.matrix(di.net.adj.matrix))
  {
    base::stop("Error in CompareNet di.net.adj.matrix is not a matrix")
  }
  if(!base::is.matrix(cmi.net.adj.matrix))
  {
    base::stop("Error in CompareNet cmi.net.adj.matrix is not a matrix")
  }
  for(i in 1:num.node)
  {
    for(j in 1:num.node)
    {
      if(di.net.adj.matrix[i,j] != cmi.net.adj.matrix[i,j])
      {
        return(0)
      }
    }
  }
  return (1)
}

################################################################################################3
#' Given two di network adjacency matrices, it prints the common edges in an output file
#'
#' Given two di nets' adjacency matrices (must of of same dim, same rownames, same colnames
#' where rows = source nodes, cols = tgt nodes, 0 and 1 represent absence and presence of
#' an edge, resp.), this function prints the common edges in an output text file.
#'
#' @param di.net.adj.matrix1 first di network adjacency matrix
#' @param di.net.adj.matrix2 second di network adjacency matrix
#' @param output.dirname output directory to store files
#'
#' @export
Print.common.di.edges <- function(di.net.adj.matrix1, di.net.adj.matrix2,  output.dirname="./OUTPUT")
{
  if(!base::is.matrix(di.net.adj.matrix1))
  {
    base::stop("Error in Print.common.di.edges di.net.adj.matrix1 is not a matrix")
  }
  if(!base::is.matrix(di.net.adj.matrix2))
  {
    base::stop("Error in Print.common.di.edges di.net.adj.matrix2 is not a matrix")
  }
  if ((dim(di.net.adj.matrix1)[1] != dim(di.net.adj.matrix2)[1]) |
      (dim(di.net.adj.matrix1)[2] != dim(di.net.adj.matrix2)[2]))
  {
    base::stop('Dimensions are not comparable!')
  }

  output.filename <- base::paste(output.dirname, 'common.di.edgelist.txt', sep = '/')
  output.file <- base::file(output.filename, 'w')
  base::rm(output.filename)

  base::cat('source', 'target', '\n', file = output.file, sep = '\t')

  edge.cnt <- 0
  for (row.idx in 1:nrow(di.net.adj.matrix1))
  {
    for (col.idx in 1:ncol(di.net.adj.matrix1))
    {
      if ((di.net.adj.matrix1[row.idx, col.idx] == 1) & (di.net.adj.matrix2[row.idx, col.idx] == 1))
      {
        base::cat(base::rownames(di.net.adj.matrix1)[row.idx], base::colnames(di.net.adj.matrix1)[col.idx], '\n',
            file = output.file, sep = '\t')

        edge.cnt <- edge.cnt + 1
      }
    }
  }

  base::close(output.file)

  base::print(edge.cnt)
}
