#' Compute Entropy matrix from the input data
#'
#' @import stats
#'
#' @param input.data input data matrix
#'
#' @return entropy matrix
#' @examples
#' df = data.frame(c(2,3,5), c(1,3,2), c(2,31,4))
#' ComputEntropy(df)
#' @export
ComputEntropy <- function(input.data)
{
  if(!base::is.data.frame(input.data))
  {
    base::stop("Error in ComputEntropy. input.data is not a data frame")
  }
  n <- base::ncol(input.data)
  entropy.matrix <- base::matrix(0, nrow = 1, ncol = n)

  for(i in 1:n)
  {
     c1 <-  stats::var(input.data[,i])
     entropy.matrix[i] <- .5*base::log((2*pi*base::exp(1))*c1)
  }
  return(entropy.matrix)
}

############################################################################################

#' Learn the mi network structure
#'
#' @import stats
#'
#' @param mut.info.matrix matrix containing mut info
#' @param mi.net.adj.matrix mi network adjacency matrix
#' @param entropy.matrix matrix containing the entropy information
#' @param alpha parameter value alpha
#'
#' @return mi network adjacency matrix
#'
#' @export
LearnMiNetStructZstat <- function(mut.info.matrix, mi.net.adj.matrix, entropy.matrix, alpha)
{
  if(!base::is.matrix(mut.info.matrix))
  {
    base::stop("Error in LearnMiNetStructZstat mut.info.matrix is not a matrix")
  }
  if(!base::is.matrix(mi.net.adj.matrix))
  {
    base::stop("Error in LearnMiNetStructZstat mi.net.adj.matrix is not a matrix")
  }
  if(!base::is.matrix(entropy.matrix))
  {
    base::stop("Error in LearnMiNetStructZstat entropy.matrix is not a matrix")
  }

  n <- base::ncol(mut.info.matrix)

  threshold <- stats::qnorm(1-(alpha/2))

  for(i in 1:(n-1))
  {
    for(j in (i+1):n)
    {
      value <- mut.info.matrix[i,j]/(entropy.matrix[i]+entropy.matrix[j])
       fisher_transform <- .5*base::log((1+value)/(1-value))
        value1 <- base::sqrt(n - 3) * base::abs(fisher_transform)

        if(value1 > threshold)
        {
          mi.net.adj.matrix[i,j] <- 1
          mi.net.adj.matrix[j,i] <- 1
        }
    }
  }
  return(mi.net.adj.matrix)
}

############################################################################################

#' Learn the mi network structure
#'
#' @import stats
#'
#' @param mut.info.matrix matrix containing mut info
#' @param mi.net.adj.matrix mi network adjacency matrix
#' @param num.nodes number of nodes
#'
#' @return mi network adjacency matrix
#' @examples
#' LearnMiNetStructRowMedian(
#' + matrix(c(0.1,0.5,0.53,0.76,0,0.12,0.43,0.65,0.23),nrow=3),
#' + matrix(c(1,0,1,0,0,0,1,0,1),nrow=3),
#' + 3)
#'
#' @export
LearnMiNetStructRowMedian <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes)
{
  if(!base::is.matrix(mut.info.matrix))
  {
    base::stop("Error in LearnMiNetStructRowMedian mut.info.matrix is not a matrix")
  }
  if(!base::is.matrix(mi.net.adj.matrix))
  {
    base::stop("Error in LearnMiNetStructRowMedian mi.net.adj.matrix is not a matrix")
  }
  for (rowIdx in 1:num.nodes)
  {
    threshold <- stats::median(mut.info.matrix[rowIdx, -rowIdx])

    for (colIdx in 1:num.nodes)
    {
      if ((colIdx != rowIdx) & (mut.info.matrix[rowIdx, colIdx] >= threshold))
      {
        mi.net.adj.matrix[rowIdx, colIdx] <- 1
      }
    }
  }

  return(mi.net.adj.matrix)
}

############################################################################################
#' Learns the CLR network Replaces all non-zero edge weights with 1.
#'
#' @import minet
#'
#' @param mut.info.matrix matrix containing mut info
#' @param mi.net.adj.matrix mi network adjacency matrix
#' @param num.nodes number of nodes
#' @param output.dirname output directory to store files
#'
#' @return mi network adjacency matrix
#'
#' @export
LearnMiNetStructClr <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes, output.dirname="./OUTPUT")
{
  if(!base::is.matrix(mut.info.matrix))
  {
    base::stop("Error in LearnMiNetStructClr mut.info.matrix is not a matrix")
  }
  if(!base::is.matrix(mi.net.adj.matrix))
  {
    base::stop("Error in LearnMiNetStructClr mi.net.adj.matrix is not a matrix")
  }
  mi.net.adj.matrix.wt <- minet::clr(mut.info.matrix) # weighted adj matrix

  # Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
  mi.net.adj.matrix.wt[base::is.nan(mi.net.adj.matrix.wt)] <- 0

  # writeLines('\n mi.net.adj.matrix.wt = \n')
  # print(mi.net.adj.matrix.wt)
  base::save(mi.net.adj.matrix.wt, file = base::paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))

  for (rowIdx in 1:num.nodes)
  {
    for (colIdx in 1:num.nodes)
    {
      if (mi.net.adj.matrix.wt[rowIdx, colIdx] != 0)
      {
        mi.net.adj.matrix[rowIdx, colIdx] <- 1
      }
    }
  }

  return(mi.net.adj.matrix)
}

############################################################################################
#' Learns CLR network from a given discretized dataset.
#'
#' Learns CLR net from a given discretized dataset. For each node, retains top 'max.fanin' number of neighbours w.r.t. edge weight
#' and removes rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
#' If there are less than that number of edges for a node, then retain all its neighbours.
#'
#' @import minet
#'
#' @param input.data.discr Discretized dataset
#' @param num.nodes number of nodes
#' @param node.names name of nodes
#' @param num.timepts number of timepoints
#' @param max.fanin number of top neighbours to retain
#' @param output.dirname output directory to store files
#'
#' @return mi network adjacency matrix
#'
#' @export
LearnClrNetFromDiscrData <- function(input.data.discr, num.nodes, node.names, num.timepts,max.fanin, output.dirname = "./OUTPUT")
{
  if(!base::is.data.frame(input.data.discr))
  {
    base::stop("Error in LearnClrNetFromDiscrData input.data.discr is not a data frame")
  }
  # Initialize mutual information matrix with zeroes
  mut.info.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = base::c(base::list(node.names), base::list(node.names)))

  num.time.series <- (base::nrow(input.data.discr) / num.timepts)

  ## Build mutual information matrix
  for (col.idx in 1:(num.nodes - 1)) {
    for (col.idx.2 in (col.idx + 1):num.nodes) {
      mut.info <- computeCmi(input.data.discr[, col.idx], input.data.discr[, col.idx.2])
      mut.info.matrix[col.idx, col.idx.2] <- mut.info
      mut.info.matrix[col.idx.2, col.idx] <- mut.info
    }
    base::rm(col.idx.2)
  }
  base::rm(col.idx)

  ## Estimate weighted adj matrix of the CLR net
  ## from the mutual info matrix
  mi.net.adj.matrix.wt <- minet::clr(mut.info.matrix)

    ## Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
    mi.net.adj.matrix.wt.curr.series[is.nan(mi.net.adj.matrix.wt.curr.series)] <- 0

    ## 'mi.net.adj.matrix.wt.curr.series' has non-neg values.
    ## Summing up all time-series-specific 'mi.net.adj.matrix.wt.curr.series' matrices
    ## produces a ('num.nodes' by 'num.nodes') matrix of non-neg
    ## values, namely 'mi.net.adj.matrix.wt'. Later 'mi.net.adj.matrix.wt'
    ## will be divided by 'num.time.series', thus producing the
    ## arthmetic mean of all time-series-specific
    ## 'mi.net.adj.matrix.wt.curr.series' matrices.
    mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt + mi.net.adj.matrix.wt.curr.series)

  ## Arthmetic mean of all time-series-specific
  ## 'mi.net.adj.matrix.wt' matrices.
  mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt / num.time.series)

  base::save(mi.net.adj.matrix.wt, file = base::paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))

  ##############################################################
  ## Begin:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################

  ## For each target node
  for (col.idx in 1:num.nodes) {

    ## Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]

    ## Count number of neighbours having positive edge weight
    num.nbrs <- base::length(edge.wts[edge.wts > 0])

    if (num.nbrs >= max.fanin) {

      ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      ## Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- base::sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]

      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1

    } else if (num.nbrs < max.fanin) {

      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  base::rm(col.idx)

  ##############################################################
  ## End:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################

  return(mi.net.adj.matrix)

}

############################################################################################
#' Learns CLR network
#'
#' Learns CLR net. For each node, retains top 'max.fanin' number of neighbours w.r.t. edge weight
#' and removes rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
#' If there are less than that number of edges for a node, then retain all its neighbours.
#'
#' @import minet
#'
#' @param mut.info.matrix matrix containing mut info
#' @param mi.net.adj.matrix mi network adjacency matrix
#' @param num.nodes number of nodes
#' @param max.fanin number of top neighbours to retain
#' @param output.dirname output directory to store files
#'
#' @return mi network adjacency matrix
#'
#' @export
LearnClrNetMfi <- function(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin, output.dirname = "./OUTPUT")
{
  if(!base::is.matrix(mut.info.matrix))
  {
    base::stop("Error in LearnClrNetMfi mut.info.matrix is not a matrix")
  }
  if(!base::is.matrix(mi.net.adj.matrix))
  {
    base::stop("Error in LearnClrNetMfi mi.net.adj.matrix is not a matrix")
  }
  mi.net.adj.matrix.wt <- minet::clr(mut.info.matrix) # weighted adj matrix

  # Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
  mi.net.adj.matrix.wt[is.nan(mi.net.adj.matrix.wt)] <- 0

  # writeLines('\n mi.net.adj.matrix.wt = \n')
  # print(mi.net.adj.matrix.wt)
  base::save(mi.net.adj.matrix.wt, file = base::paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))

  # For each target node
  for (col.idx in 1:num.nodes)
  {
    # Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]

    # Count number of neighbours having positive edge weight
    num.nbrs <- base::length(edge.wts[edge.wts > 0])

    if (num.nbrs >= max.fanin)
    {
      # Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      # Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- base::sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]

      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1

      ## The following line is not required since 'mi.net.adj.matrix' is initialized
      ## with all zeroes
      # mi.net.adj.matrix[-(valid.nbrs), col.idx] <- 0
    }
    else if (num.nbrs < max.fanin)
    {
      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }

  return(mi.net.adj.matrix)
}

############################################################################################
#' Learns CLR2 network
#'
#' Learns CLR2 net. For each node, retains top 'max.fanin' number of neighbours w.r.t. edge weight
#' and removes rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
#' If there are less than that number of edges for a node, then retain all its neighbours.
#'
#' @import minet
#'
#' @param input.data.discr Discretized dataset
#' @param num.nodes number of nodes
#' @param node.names name of nodes
#' @param num.timepts number of timepoints
#' @param max.fanin number of top neighbours to retain
#' @param output.dirname output directory name
#' @param mi.net.adj.matrix mi network adjacency matrix
#'
#' @return mi network adjacency matrix
#'
#' @export
LearnClr2NetMfi <- function(input.data.discr, num.nodes, node.names, num.timepts,
                            max.fanin, output.dirname = "./OUTPUT", mi.net.adj.matrix)
{
  if(!base::is.data.frame(input.data.discr))
  {
    base::stop("Error in LearnClr2NetMfi input.data.discr is not a data frame")
  }
  if(!base::is.matrix(mi.net.adj.matrix))
  {
    base::stop("Error in LearnClr2NetMfi mi.net.adj.matrix is not a matrix")
  }
  ## Initialize weighted adjacency matrix of the mutual information network
  mi.net.adj.matrix.wt <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                              dimnames = base::c(base::list(node.names), base::list(node.names)))

  ## Total number of time series
  ## = number of measurements (replicates) per time pt
  ## = (total number of measurements / number of time pts).
  num.time.series <- (base::nrow(input.data.discr) / num.timepts)

  for (time.series.idx in 1:num.time.series) {

    ## First time point of the current time series
    first.time.pt.curr.series <- (((time.series.idx - 1) * num.timepts) + 1)

    ## Last time point of the current time series
    last.time.pt.curr.series <- (time.series.idx * num.timepts)

    ## Discretized data of the current time series
    input.data.discr.curr.series <- input.data.discr[first.time.pt.curr.series:last.time.pt.curr.series, ]
    base::rm(first.time.pt.curr.series, last.time.pt.curr.series)

    # Initialize mutual information matrix with zeroes
    mut.info.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = base::c(base::list(node.names), base::list(node.names)))

    ## Build mutual information matrix
    for (col.idx in 1:(num.nodes - 1)) {
      for (col.idx.2 in (col.idx + 1):num.nodes) {

        ## compute_cmi.R
        mut.info <- computeCmi(input.data.discr.curr.series[, col.idx], input.data.discr.curr.series[, col.idx.2])

        mut.info.matrix[col.idx, col.idx.2] <- mut.info
        mut.info.matrix[col.idx.2, col.idx] <- mut.info
      }
      base::rm(col.idx.2)
    }
    base::rm(col.idx)

    ## Estimate weighted adj matrix of the CLR net
    ## corr. to the current time series
    mi.net.adj.matrix.wt.curr.series <- minet::clr(mut.info.matrix)

    ## Replace 'NaN' with zero. 'NaN' is produced when a corr. variable has variance zero.
    mi.net.adj.matrix.wt.curr.series[is.nan(mi.net.adj.matrix.wt.curr.series)] <- 0

    ## 'mi.net.adj.matrix.wt.curr.series' has non-neg values.
    ## Summing up all time-series-specific 'mi.net.adj.matrix.wt.curr.series' matrices
    ## produces a ('num.nodes' by 'num.nodes') matrix of non-neg
    ## values, namely 'mi.net.adj.matrix.wt'. Later 'mi.net.adj.matrix.wt'
    ## will be divided by 'num.time.series', thus producing the
    ## arthmetic mean of all time-series-specific
    ## 'mi.net.adj.matrix.wt.curr.series' matrices.
    mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt + mi.net.adj.matrix.wt.curr.series)

  }
  base::rm(time.series.idx)

  ## Arthmetic mean of all time-series-specific
  ## 'mi.net.adj.matrix.wt' matrices.
  mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt / num.time.series)

  base::save(mi.net.adj.matrix.wt, file = base::paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))

  ##############################################################
  ## Begin:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################

  ## For each target node
  for (col.idx in 1:num.nodes) {

    ## Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]

    ## Count number of neighbours having positive edge weight
    num.nbrs <- base::length(edge.wts[edge.wts > 0])

    if (num.nbrs >= max.fanin) {

      ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      ## Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- base::sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]

      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1

    } else if (num.nbrs < max.fanin) {

      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  base::rm(col.idx)

  ##############################################################
  ## End:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################

  return(mi.net.adj.matrix)

}

############################################################################################

############################################################################################
#' Learn CLR2.1 network
#'
#' Learns CLR2.1 net. For each node, retains top 'max.fanin' number of neighbours w.r.t. edge weight
#' and removes rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
#' If there are less than that number of edges for a node, then retain all its neighbours.
#'
#' @import stats
#'
#' @param input.data.discr Discretized dataset
#' @param num.nodes number of nodes
#' @param node.names name of nodes
#' @param num.timepts number of timepoints
#' @param max.fanin number of top neighbours to retain
#' @param output.dirname output directory name
#' @param mi.net.adj.matrix mi network adjacency matrix
#'
#' @return mi network adjacency matrix
#'
#' @export
LearnClrNetMfiVer2.1 <- function(input.data.discr, num.nodes, node.names, num.timepts,
                                 max.fanin, output.dirname = "./OUTPUT", mi.net.adj.matrix)
{
  if(!base::is.data.frame(input.data.discr))
  {
    base::stop("Error in LearnClrNetMfiVer2.1 input.data.discr is not a data frame")
  }
  if(!base::is.matrix(mi.net.adj.matrix))
  {
    base::stop("Error in LearnClrNetMfiVer2.1 mi.net.adj.matrix is not a matrix")
  }
  ## Initialize weighted adjacency matrix of the mutual information network
  mi.net.adj.matrix.wt <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                                 dimnames = base::c(base::list(node.names), base::list(node.names)))

  ## Total number of time series
  ## = number of measurements (replicates) per time pt
  ## = (total number of measurements / number of time pts).
  num.time.series <- (base::nrow(input.data.discr) / num.timepts)

  for (time.series.idx in 1:num.time.series) {

    ## First time point of the current time series
    first.time.pt.curr.series <- (((time.series.idx - 1) * num.timepts) + 1)

    ## Last time point of the current time series
    last.time.pt.curr.series <- (time.series.idx * num.timepts)

    ## Discretized data of the current time series
    input.data.discr.curr.series <- input.data.discr[first.time.pt.curr.series:last.time.pt.curr.series, ]
    base::rm(first.time.pt.curr.series, last.time.pt.curr.series)

    # Initialize mutual information matrix with zeroes
    mut.info.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes, dimnames = c(list(node.names), list(node.names)))

    ## Build assymetric mutual information matrix
    for (row.idx in 1:num.nodes) {
      for (col.idx in 1:num.nodes) {

        ## compute_cmi.R
        mut.info <- computeCmi(input.data.discr.curr.series[1:(num.timepts-1), row.idx],
                               input.data.discr.curr.series[2:num.timepts, col.idx])

        mut.info.matrix[row.idx, col.idx] <- mut.info
      }
      base::rm(col.idx)
    }
    base::rm(row.idx)

    # print('mut.info.matrix')
    # print(time.series.idx)
    # print(mut.info.matrix)

    ## Initialize weighted adj matrix of the CLR net
    ## corr. to the current time series
    mi.net.adj.matrix.wt.curr.series <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                                               dimnames = c(list(node.names), list(node.names)))
    ## Compute 'mi.net.adj.matrix.wt.curr.series[src.node.idx, tgt.node.idx]'
    for (src.node.idx in 1:num.nodes) {
      for (tgt.node.idx in 1:num.nodes) {

        if (mut.info.matrix[src.node.idx, tgt.node.idx] == .Machine$double.xmax) {
          ## Perfectly correlated

          mi.net.adj.matrix.wt.curr.series[src.node.idx, tgt.node.idx] <- .Machine$double.xmax

        } else {

          src.clr.mean <- base::mean(mut.info.matrix[src.node.idx, ])
          # if (src.clr.mean > .Machine$double.xmax) {
          #   src.clr.mean <- .Machine$double.xmax
          # }

          src.clr.sd <- stats::sd(mut.info.matrix[src.node.idx, ])
          # if (src.clr.sd > .Machine$double.xmax) {
          #   src.clr.sd <- .Machine$double.xmax
          # }

          tgt.clr.mean <- base::mean(mut.info.matrix[, tgt.node.idx])
          # if (tgt.clr.mean > .Machine$double.xmax) {
          #   tgt.clr.mean <- .Machine$double.xmax
          # }

          tgt.clr.sd <- stats::sd(mut.info.matrix[, tgt.node.idx])
          # if (tgt.clr.sd > .Machine$double.xmax) {
          #   tgt.clr.sd <- .Machine$double.xmax
          # }

          if ((src.clr.mean >= .Machine$double.xmax) | (src.clr.sd >= .Machine$double.xmax) |
              (tgt.clr.mean >= .Machine$double.xmax) | (tgt.clr.sd >= .Machine$double.xmax)) {

            ## There exists far better candidate source node(s)
            mi.net.adj.matrix.wt.curr.series[src.node.idx, tgt.node.idx] <- 0

          } else {

            tmp <- 0
            if (src.clr.sd != 0) {
              tmp <- (mut.info.matrix[src.node.idx, tgt.node.idx] - src.clr.mean) / src.clr.sd
            }
            z.src <- base::max(0, tmp)

            ## Re-initilize in case 'tgt.clr.sd == 0'
            tmp <- 0
            if (tgt.clr.sd != 0) {
              tmp <- (mut.info.matrix[src.node.idx, tgt.node.idx] - tgt.clr.mean) / tgt.clr.sd
            }
            z.tgt <- base::max(0, tmp)

            base::rm(tmp)

            clr.edge.wt <- ((z.tgt)^2 + (z.src)^2)
            clr.edge.wt <- base::sqrt(clr.edge.wt)

            mi.net.adj.matrix.wt.curr.series[src.node.idx, tgt.node.idx] <- clr.edge.wt
          }
        }

      }
      base::rm(tgt.node.idx)
    }
    base::rm(src.node.idx)

    ## 'mi.net.adj.matrix.wt.curr.series' has non-neg values.
    ## Summing up all time-series-specific 'mi.net.adj.matrix.wt.curr.series' matrices
    ## produces a ('num.nodes' by 'num.nodes') matrix of non-neg
    ## values, namely 'mi.net.adj.matrix.wt'. Later 'mi.net.adj.matrix.wt'
    ## will be divided by 'num.time.series', thus producing the
    ## arthmetic mean of all time-series-specific
    ## 'mi.net.adj.matrix.wt.curr.series' matrices.
    mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt + mi.net.adj.matrix.wt.curr.series)

    ## '.Machine$double.xmax + .Machine$double.xmax = Inf'
    ## Prevent 'Inf'
    mi.net.adj.matrix.wt[mi.net.adj.matrix.wt > .Machine$double.xmax] <- .Machine$double.xmax

  }
  base::rm(time.series.idx)

  ## Arthmetic mean of all time-series-specific
  ## 'mi.net.adj.matrix.wt' matrices.
  mi.net.adj.matrix.wt <- (mi.net.adj.matrix.wt / num.time.series)

  base::save(mi.net.adj.matrix.wt, file = base::paste(output.dirname, 'mi.net.adj.matrix.wt.RData', sep = '/'))

  ##############################################################
  ## Begin:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################

  ## For each target node
  for (col.idx in 1:num.nodes) {

    ## Weights of the edges with the target node
    edge.wts <- mi.net.adj.matrix.wt[, col.idx]

    ## Count number of neighbours having positive edge weight
    num.nbrs <- base::length(edge.wts[edge.wts > 0])

    if (num.nbrs >= max.fanin) {

      ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
      ## Tie is broken in favour of the neighbour having smaller index.
      valid.nbrs <- base::sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]

      mi.net.adj.matrix[valid.nbrs, col.idx] <- 1

    } else if (num.nbrs < max.fanin) {

      # Retain all the neighbours
      mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
    }
  }
  base::rm(col.idx)

  ##############################################################
  ## End:
  ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
  ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
  ## and 'max.fanin'
  ##############################################################

  return(mi.net.adj.matrix)

}

############################################################################################

############################################################################################
#' Learn CLR3 network
#'
#' Learns CLR3 net. For each node, retains top 'max.fanin' number of neighbours w.r.t. edge weight
#' and removes rest of the edges. Tie is broken in favour of the neighbour having smaller node index.
#' If there are less than that number of edges for a node, then retain all its neighbours.
#'
#' @import stats
#'
#' @param input.data.discr.3D 3D Discretized dataset
#' @param num.nodes number of nodes
#' @param node.names name of nodes
#' @param num.timepts number of timepoints
#' @param max.fanin number of top neighbours to retain
#' @param mi.net.adj.matrix.list mi network adjacency matrix list
#'
#' @return mi network adjacency matrix list
#'
#' @export
LearnClr3NetMfi <- function(input.data.discr.3D, num.nodes, node.names, num.timepts,
                            max.fanin, mi.net.adj.matrix.list)
{
  if(!base::is.array(input.data.discr.3D))
  {
    base::stop("Error in LearnClr3NetMfi input.data.discr.3D is not an array")
  }
  if(!base::is.list(mi.net.adj.matrix.list))
  {
    base::stop("Error in LearnClr3NetMfi mi.net.adj.matrix.list is not a matrix")
  }
  ## Here, each 'time.pt.idx' represents time interval
  ## ('time.pt.idx', ('time.pt.idx' + 1))
  for (time.pt.idx in 1:(num.timepts - 1)) {

    ## Discretized data corr. to the current time interval
    input.data.discr.3D.curr.ival <- input.data.discr.3D[(time.pt.idx:(time.pt.idx + 1)), , ]

    candidate.parent.node.names <- base::c()
    candidate.tgt.node.names <- base::c()
    for (curr.node.name in node.names) {
      parent.full.name <- base::paste(curr.node.name, as.character(time.pt.idx), sep = '_t')
      candidate.parent.node.names <- base::c(candidate.parent.node.names, parent.full.name)
      base::rm(parent.full.name)

      tgt.full.name <- base::paste(curr.node.name, as.character(time.pt.idx + 1), sep = '_t')
      candidate.tgt.node.names <- base::c(candidate.tgt.node.names, tgt.full.name)
      base::rm(tgt.full.name)
    }
    base::rm(curr.node.name)

    ## Initialize mutual information matrix with zeroes
    mut.info.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                              dimnames = c(list(candidate.parent.node.names),
                                           list(candidate.tgt.node.names)))

    ## Build mutual information matrix
    for (col.idx in 1:(num.nodes - 1)) {
      for (col.idx.2 in (col.idx + 1):num.nodes) {

        ## compute_cmi.R
        ## (dim1 == 1) => time.pt.idx
        ## (dim1 == 2) => (time.pt.idx + 1)
        mut.info <- computeCmi(input.data.discr.3D.curr.ival[1, col.idx, ],
                               input.data.discr.3D.curr.ival[2, col.idx.2, ])

        mut.info.matrix[col.idx, col.idx.2] <- mut.info
        mut.info.matrix[col.idx.2, col.idx] <- mut.info
      }
      base::rm(col.idx.2)
    }
    base::rm(col.idx)

    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix.wt <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                              dimnames = base::c(base::list(candidate.parent.node.names),
                                                 base::list(candidate.tgt.node.names)))

    candidate.parent.mean.sd <- base::matrix(0, nrow = num.nodes, ncol = 2)
    base::rownames(candidate.parent.mean.sd) <- candidate.parent.node.names
    base::colnames(candidate.parent.mean.sd) <- base::c('clr.mean', 'clr.sd')

    ## Calculate sample mean and sample standard deviation of the given nodes.
    ## It is calculated acc. to the logic in function 'clr()' in
    ## R package 'minet' (version 3.36.0).
    ## The aforementioned 'clr()' function uses 'clr.cpp' to perform
    ## the calculation. Here, the same logic is re-implemented in R.
    ## Since the in-built 'mean()' and 'sd()' functions in
    ## R version 3.3.2 follows the exact same logic, therefore, the
    ## re-implementation is straight-forward.
    for (parent.name in candidate.parent.node.names) {
      ## arithmetic mean
      candidate.parent.mean.sd[parent.name, 'clr.mean'] <- base::mean(mut.info.matrix[parent.name, ])

      ## var <- 0
      ## for (each sample) {
      ##  sd <- (mean - sample val)
      ##  var <- var + (sd^2)
      ## }
      ## var <- var / (n -1) ## where n = number of samples
      ## sd <- sqrt(var)
      candidate.parent.mean.sd[parent.name, 'clr.sd'] <- stats::sd(mut.info.matrix[parent.name, ])
    }
    base::rm(parent.name)

    for (tgt.node.name in candidate.tgt.node.names) {

      tgt.clr.mean <- base::mean(mut.info.matrix[, tgt.node.name])
      tgt.clr.sd <- stats::sd(mut.info.matrix[, tgt.node.name])

      ## Edge weights of the CLR net are calculated acc. to
      ## 'minet::clr()'
      for (candidate.parent.name in candidate.parent.node.names) {

        tmp <- 0
        if (tgt.clr.sd != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] - tgt.clr.mean) / tgt.clr.sd
        }
        z.tgt <- base::max(0, tmp)

        tmp <- 0
        if (candidate.parent.mean.sd[candidate.parent.name, 'clr.sd'] != 0) {
          tmp <- (mut.info.matrix[candidate.parent.name, tgt.node.name] -
                    candidate.parent.mean.sd[candidate.parent.name, 'clr.mean']) /
            candidate.parent.mean.sd[candidate.parent.name, 'clr.sd']
        }
        z.parent <- base::max(0, tmp)

        base::rm(tmp)

        clr.edge.wt <- ((z.tgt)^2 + (z.parent)^2)
        clr.edge.wt <- base::sqrt(clr.edge.wt)
        base::rm(z.tgt, z.parent)

        mi.net.adj.matrix.wt[candidate.parent.name, tgt.node.name] <- clr.edge.wt
      }
      base::rm(candidate.parent.name)

    }
    base::rm(tgt.node.name)

    ##############################################################
    ## Begin:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################

    ## Initialize unweighted adjacency matrix of the CLR net
    ## corr. to the current time interval
    mi.net.adj.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                                   dimnames = base::c(base::list(candidate.parent.node.names),
                                                      base::list(candidate.tgt.node.names)))

    ## For each target node
    for (col.idx in 1:num.nodes) {

      ## Weights of the edges with the target node
      edge.wts <- mi.net.adj.matrix.wt[, col.idx]

      ## Count number of neighbours having positive edge weight
      num.nbrs <- base::length(edge.wts[edge.wts > 0])

      if (num.nbrs >= max.fanin) {

        ## Return indices of the top 'max.fanin' number of neighbours w.r.t. edge weight.
        ## Tie is broken in favour of the neighbour having smaller index.
        valid.nbrs <- base::sort(edge.wts, decreasing = TRUE, index.return = TRUE)$ix[1:max.fanin]

        mi.net.adj.matrix[valid.nbrs, col.idx] <- 1

      } else if (num.nbrs < max.fanin) {

        # Retain all the neighbours
        mi.net.adj.matrix[edge.wts > 0, col.idx] <- 1
      }
    }
    base::rm(col.idx)

    ##############################################################
    ## End:
    ## Estimate the unweighted adjacency matrix 'mi.net.adj.matrix'
    ## using the weighted adjacency matrix 'mi.net.adj.matrix.wt'
    ## and 'max.fanin'
    ##############################################################


    mi.net.adj.matrix.list[[time.pt.idx]] <- mi.net.adj.matrix
  }
  base::rm(time.pt.idx)

  return(mi.net.adj.matrix.list)
}

############################################################################################
