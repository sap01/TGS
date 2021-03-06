#' Given a directed network, convert it into an undirected network
#'
#' @param di.net a directed network
#'
#' @return directed network converted to undirected network
# #' @examples
#' ConvertDinetToUndinet(matrix(c(0,1,0,0),nrow=2))
#'
#' @keywords internal
#' @noRd
ConvertDinetToUndinet <- function(di.net) {
  if(!base::is.matrix(di.net))
  {
    base::stop("Error in ConvertDinetToUndinet di.net is not a matrix")
  }
  ## Initilize undirected network with zeroes
  undi.net <- di.net
  undi.net[undi.net] <- 0

  for (row.idx in 1:nrow(di.net)) {
    for (col.idx in 1:ncol(di.net)) {
      if (di.net[row.idx, col.idx] == 1) {
        undi.net[row.idx, col.idx] <- 1
        undi.net[col.idx, row.idx] <- 1
      }
    }
    base::rm(col.idx)
  }
  base::rm(row.idx)

return(undi.net)
}
