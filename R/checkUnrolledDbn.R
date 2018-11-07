#' Checks whether the given unrolled DBN follows 1st Markov order or not
#'
#' @param unrolled.DBN.adj.matrix An unrolled DBN adjacency matrix
#'
#' @return whether the given unrolled DBN follows 1st Markov order or not
#'
#' @examples
#' checkUnrolledDbn(matrix(c(0,0,0,0),nrow=2,ncol=2))
#'
#' @export
checkUnrolledDbn <- function(unrolled.DBN.adj.matrix)
{
  if(!base::is.matrix(unrolled.DBN.adj.matrix))
  {
    base::stop("Error in checkUnrolledDbn. unrolled.DBN.adj.matrix is not a matrix")
  }
  for(row.idx in 1:nrow(unrolled.DBN.adj.matrix))
  {
    for(col.idx in 1:ncol(unrolled.DBN.adj.matrix))
    {
      if (unrolled.DBN.adj.matrix[row.idx, col.idx] == 1)
      {
        src.node.name <- base::rownames(unrolled.DBN.adj.matrix)[row.idx]
        tgt.node.name <- base::colnames(unrolled.DBN.adj.matrix)[col.idx]

        src.time.pt.name <- base::strsplit(src.node.name, '_', fixed = TRUE)[[1]][2]
        src.time.pt <- base::substr(src.time.pt.name, 2, nchar(src.time.pt.name))
        src.time.pt <- base::as.integer(src.time.pt)

        tgt.time.pt.name <- base::strsplit(tgt.node.name, '_', fixed = TRUE)[[1]][2]
        tgt.time.pt <- base::substr(tgt.time.pt.name, 2, nchar(tgt.time.pt.name))
        tgt.time.pt <- base::as.integer(tgt.time.pt)

        if (src.time.pt != (tgt.time.pt - 1))
        {
          base::cat(src.node.name,'\t',tgt.node.name,'\n')
        }
      }
    }
  }
}
