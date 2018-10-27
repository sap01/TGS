#' @title Convert a given unrolled Dynamic Bayesian Network (DBN) into a rolled DBN using different rolling methods
#'
#' Rolls time-varying networks into a single time-invariant network.
#' This function is compatible with the time-varying networks learnt through
#' learnDbnStruct3dParDeg1.R::learnDbnStructMo1Layer3dParDeg1().
#'
#' @param num.nodes Number of the desired nodes in the rolled DBN
#' @param node.names Names of the desired nodes in the rolled DBN
#' @param num.timepts Number of time points in the unrolled DBN
#' @param unrolled.DBN.adj.matrix Given unrolled DBN adjacency matrix. It is a 2D matrix of dimension ((num.nodes X num.timepts) X (num.nodes X num.timepts)).
#' @param roll.method Which rolling method to use from {'any', 'all', or some real number in (0, 1), like - 0.5}.
#' @param allow.self.loop Boolean to decide whether to allow self loop or not in the rolled DBN
#'
#' @return rolled.DBN.adj.matrix Return the rolled DBN adjacency matrix. It is a 2D matrix of dimension (num.nodes * num.nodes).
#'
#' @export
rollDbn <- function(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, roll.method, allow.self.loop)
{
  if(!base::is.matrix(unrolled.DBN.adj.matrix))
  {
    base::stop("Error in rollDbn unrolled.DBN.adj.matrix is not a matrix")
  }
  ## Initialize rolled DBN adj matrix as a zero matrix
  rolled.DBN.adj.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                                  dimnames = base::c(base::list(node.names), base::list(node.names)))

  num.time.trans <- num.timepts - 1 # Num of time transitions

  # todo: replace with foreach dopar
  for (tgt.node.idx in 1:ncol(rolled.DBN.adj.matrix))
  {
    tgt.node.name <- base::colnames(rolled.DBN.adj.matrix)[tgt.node.idx]

    # grep('^G1', tmpvec, fixed = FALSE) returns the indices of the elements in 'tmpvec' whose values start with 'G1'.
    # '^G1' is the given pattern.
    # 'fixed = FALSE' represents that the given pattern is a regular expression.
    unrolled.DBN.tgt.node.indices <- base::grep(base::paste('^', tgt.node.name, sep = ''),
                                                base::colnames(unrolled.DBN.adj.matrix),
                                        fixed = FALSE)
    unrolled.DBN.adj.matrix.tgt.node <- unrolled.DBN.adj.matrix[, unrolled.DBN.tgt.node.indices]

    # If the value corr. to a row in 'unrolled.DBN.adj.matrix.tgt.node.single.col'
    # is greater than zero, then the node corr. to the row name is a parent of the target node
    unrolled.DBN.adj.matrix.tgt.node.single.col <-
      base::matrix(rowSums(unrolled.DBN.adj.matrix.tgt.node),
             nrow = base::nrow(unrolled.DBN.adj.matrix.tgt.node), ncol = 1,
             dimnames = base::c(base::list(base::rownames(unrolled.DBN.adj.matrix.tgt.node)),
                          tgt.node.name))

    # After execution of this for loop,
    # rolled.DBN.adj.matrix[src.node.name, tgt.node.name] represents how many times there is an edge from
    # the src node to the tgt node in the unrolled DBN. The value is an integer in the interval [0, num.time.trans].
     for (src.node.name in node.names)
    {
      unrolled.DBN.src.node.indices <- base::grep(paste('^', src.node.name, sep = ''),
                                                  base::rownames(unrolled.DBN.adj.matrix.tgt.node.single.col),
                                          fixed = FALSE)

      rolled.DBN.adj.matrix[src.node.name, tgt.node.name] <-
        base::sum(unrolled.DBN.adj.matrix.tgt.node.single.col[unrolled.DBN.src.node.indices, tgt.node.name])

     }
  }

  roll.threshold <- NULL

  if(base::is.character(roll.method))
  {
    if (roll.method == 'any') # Insert an edge in rolled DBN if it is present at least for one time transition
    {
      roll.threshold <- 1

    }
    else if (roll.method == 'all') # Insert an edge in rolled DBN if it is present at every time transition
    {
      roll.threshold <- num.time.trans
    }
  }
  else if (base::is.numeric(roll.method)) # Insert an edge in rolled DBN if it is present at at least (roll.method * num.time.trans) number of time transitions
  {
    if ((roll.method > 0) & (roll.method < 1))
    {
      roll.threshold <- num.time.trans * roll.method
    }
    else
    {
      # print('\'roll.method\' accepts numeric values in the interval (0,1)')
      base::stop('\'roll.method\' accepts numeric values in the interval (0,1)')
    }
  }

  # writeLines('\n rolled.DBN.adj.matrix = \n')
  # print(rolled.DBN.adj.matrix)

  for (tgt.node.idx in 1:ncol(rolled.DBN.adj.matrix))
  {
    for (src.node.idx in 1:nrow(rolled.DBN.adj.matrix))
    {
      if (rolled.DBN.adj.matrix[src.node.idx, tgt.node.idx] >= roll.threshold)
      {
        rolled.DBN.adj.matrix[src.node.idx, tgt.node.idx] <- 1
      }
      else
      {
        rolled.DBN.adj.matrix[src.node.idx, tgt.node.idx] <- 0
      }
    }
  }

  # Remove self loops if 'allow.self.loop' = FALSE
  if (!allow.self.loop)
  {
    base::diag(rolled.DBN.adj.matrix) <- 0
  }

  return(rolled.DBN.adj.matrix)
}

####################################################################################################3
#' @title Convert a given unrolled Dynamic Bayesian Network (DBN) into a rolled DBN using different rolling methods
#'
#' Rolls time-varying networks into a single time-invariant network.
#' This function is compatible with the time-varying networks learnt through
#' learnDbnStruct3dParDeg1.R::learnDbnStructMo1Layer3dParDeg1_v2().
#'
#' @param num.nodes Number of the desired nodes in the rolled DBN
#' @param node.names Names of the desired nodes in the rolled DBN
#' @param num.timepts Number of time points in the unrolled DBN
#' @param unrolled.DBN.adj.matrix.list Given time-varying network adjacency list. Its length =
#' num.time.trans = (num.timepts - 1). The t^{th} element of the list represents the predicted
#' network adjacency matrix of the t^{th} time transition. This matrix is of dimension
#' (num.nodes \ times num.nodes).
#' @param roll.method Which rolling method to use from {'any', 'all', or some real number in (0, 1), like - 0.5}.
#' @param allow.self.loop Boolean to decide whether to allow self loop or not in the rolled DBN
#'
#' @return rolled.DBN.adj.matrix Return the rolled DBN adjacency matrix. It is a 2D matrix of dimension (num.nodes * num.nodes).
#'
#' @export
rollDbn_v2 <- function(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix.list,
                       roll.method, allow.self.loop)
{
  if(!base::is.list(unrolled.DBN.adj.matrix.list))
  {
    base::stop("Error in rollDbn_v2 unrolled.DBN.adj.matrix.list is not a list")
  }
  ## Initialize rolled DBN adj matrix as a zero matrix
  rolled.DBN.adj.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                                  dimnames = base::c(base::list(node.names), base::list(node.names)))

  num.time.trans <- num.timepts - 1 # Num of time transitions

  for (list.idx in 1:num.time.trans)
  {
    rolled.DBN.adj.matrix <- rolled.DBN.adj.matrix + unrolled.DBN.adj.matrix.list[[list.idx]]
  }
  base::rm(list.idx)

  roll.threshold <- NULL

  if(base::is.character(roll.method))
  {
    if (roll.method == 'any') # Insert an edge in rolled DBN if it is present at least for one time transition
    {
      roll.threshold <- 1

    }
    else if (roll.method == 'all') # Insert an edge in rolled DBN if it is present at every time transition
    {
      roll.threshold <- num.time.trans
    }
  }
  else if (base::is.numeric(roll.method)) # Insert an edge in rolled DBN if it is present at at least (roll.method * num.time.trans) number of time transitions
  {
    if ((roll.method > 0) & (roll.method < 1))
    {
      roll.threshold <- num.time.trans * roll.method
    }
    else
    {
      # print('\'roll.method\' accepts numeric values in the interval (0,1)')
      base::stop('\'roll.method\' accepts numeric values in the interval (0,1)')
    }
  }

  # writeLines('\n rolled.DBN.adj.matrix = \n')
  # print(rolled.DBN.adj.matrix)

  for (tgt.node.idx in 1:ncol(rolled.DBN.adj.matrix))
  {
    for (src.node.idx in 1:nrow(rolled.DBN.adj.matrix))
    {
      if (rolled.DBN.adj.matrix[src.node.idx, tgt.node.idx] >= roll.threshold)
      {
        rolled.DBN.adj.matrix[src.node.idx, tgt.node.idx] <- 1
      }
      else
      {
        rolled.DBN.adj.matrix[src.node.idx, tgt.node.idx] <- 0
      }
    }
  }
  base::rm(tgt.node.idx)

  # Remove self loops if 'allow.self.loop' = FALSE
  if (!allow.self.loop)
  {
    base::diag(rolled.DBN.adj.matrix) <- 0
  }

  return(rolled.DBN.adj.matrix)
}
