#' Implementing the TGS Algorithm.
#'
#' It Uncovers the underlying temporal sequence of gene regulatory events in the form of time-varying Gene Regulatory Networks (GRNs)
#' It learns the time-varying GRN structures independently of each other, without imposing any structural constraint. However,
#' it is time intensive and hence not suitable for large-scale GRNs.
#'
#' @param isfile 1 if the parameters are given in a json file otherwise 0.
#' @param input.data.filename name of the file containing the data without the directory name. It can be .tsv or .Rdata file only.
#' @param num.timepts number of distinct timepoints
#' @param true.net.filename File containing the true network without the directory name. In case non empty then should contain .Rdata file with object name 'true.net.adj.matrix'
#' @param input.wt.data.filename File containing Input Wild Type data without the directory name. If it is not empty hen must be a .tsv file having first row containing names of genes except (1,1) and second row should have the WT values except (2,1)th cell.
#' @param is.discrete whether the data is discretized or not
#' @param num.discr.levels number of discreteized levels that each gene has (if already discretized) or it should have (if it is to be discretized)
#' @param discr.algo The algo to follow in case of not discretized. Possible values:- {discretizeData.2L.Tesla,discretizeData.2L.wt.l}
#' @param mi.estimator which method to use for estimating the mutual information matrix. Generally ‘mi.pca.cmi’ is used
#' @param apply.aracne ARACNE is applied to refine the mutual information matrix. In case true then TGS+ version is used otherwise the original TGS variant is executed
#' @param clr.algo CLR algo to use Possible values :- {CLR, CLR2, CLR2.1, CLR3, spearman}
#' @param max.fanin the maximum number of regulators each gene can have
#' @param allow.self.loop Whether to allow self loops in the graph
#' @param scoring.func Which scoring func to use
#' @param input.dirname Name of the directory where input files are. By default your current directory unless specified otherwise
#' @param output.dirname Name of the directory where output files are to be stored. By default your current directory unless specified otherwise
#' @param json.file name of the json file along with directory if parameters are to be read from json file
#'
#' @import rjson
#' @import minet
#' @import utils
#'
#' @examples
#' \dontrun{
#' LearnTgs(0,input.data.filename = "DmLc3E.RData", num.timepts = 6, is.discrete = TRUE,
#'  num.discr.levels = 2, mi.estimator = "mi.pca.cmi", apply.aracne = FALSE,
#'  clr.algo = "CLR", max.fanin = 14, allow.self.loop = TRUE,
#'  input.dirname = "location where file is stored",
#'  output.dirname = "location where output needs to be stored")
#'
#' LearnTgs(0,input.data.filename = "DmLc3L.RData", num.timepts = 2, is.discrete = TRUE,
#'  num.discr.levels = 2, mi.estimator = "mi.pca.cmi", apply.aracne = FALSE,
#'  clr.algo = "CLR", max.fanin = 14, allow.self.loop = TRUE,
#'  input.dirname = "location where file is stored",
#'  output.dirname = "location where output needs to be stored")
#' }
#'
#' @export
LearnTgs <- function(isfile = 0,
                      input.data.filename = "",
                      num.timepts = 0,
                      true.net.filename = "",
                      input.wt.data.filename = "",
                      is.discrete = TRUE,
                      num.discr.levels = 2,
                      discr.algo = "",
                      mi.estimator = "mi.pca.cmi",
                      apply.aracne = FALSE,
                      clr.algo = "CLR",
                      max.fanin = 14,
                      allow.self.loop = TRUE,
                      scoring.func = "BIC",
                      input.dirname = "",
                      output.dirname = "",
                      json.file = ""
                      )
{
  if(isfile!=0)
  {
    ##------------------------------------------------------------
    ## Begin: Read User-defined input Params
    ##------------------------------------------------------------

    input.params <- rjson::fromJSON(file = json.file)
    input.data.filename <- input.params$input.data.filename
    num.timepts <- input.params$num.timepts
    true.net.filename <- input.params$true.net.filename
    input.wt.data.filename <- input.params$input.wt.data.filename
    is.discrete <- input.params$is.discrete
    num.discr.levels <- input.params$num.discr.levels
    discr.algo <- input.params$discr.algo
    mi.estimator <- input.params$mi.estimator
    apply.aracne <- input.params$apply.aracne
    clr.algo <- input.params$clr.algo
    max.fanin <- input.params$max.fanin
    allow.self.loop <- input.params$allow.self.loop
    rm(input.params)

    ##------------------------------------------------------------
    ## End: Read User-defined input Params
    ##------------------------------------------------------------
  }

  if(input.dirname=="")
  {
    input.dirname <- base::getwd()
  }
  if(output.dirname=="")
  {
    output.dirname <- base::getwd()
  }

  input.data.filename <- base::paste(input.dirname, input.data.filename, sep = '/')
  if (true.net.filename != '')
  {
    true.net.filename <- base::paste(input.dirname, true.net.filename, sep = '/')
  }
  if (input.wt.data.filename != '')
  {
    input.wt.data.filename <- base::paste(input.dirname, input.wt.data.filename, sep = '/')
  }

  ##------------------------------------------------------------
  ## Begin: Main program
  ##------------------------------------------------------------

  base::print('The output directory name is:')
  base::print(output.dirname)
  base::print('') ## to append a blank line

  ## Save console output in a file named 'output.txt' inside the output directory.
  output.filename <- base::paste(output.dirname, 'output.txt', sep = '/')
  output.file.conn <- base::file(output.filename, open = "wt")
  base::sink(output.file.conn)

  ##------------------------------------------------------------
  ## Begin: Read input data file
  ##------------------------------------------------------------

  ## Begin: Find file extension of the input data file. Only '.tsv' and '.RData'
  ## are allowed.
  ## Split the string at every '.' and consider the last substring as the
  ## file extension.
  input.data.filename.ext <- base::unlist(base::strsplit(input.data.filename, '[.]'))
  ## End: Find file extension of the input data file. Only '.tsv' and '.RData'
  ## are allowed.

  ## Initialize input data
  input.data <- NULL
  if (input.data.filename.ext[base::length(input.data.filename.ext)] == 'tsv') {
    input.data <- utils::read.table(input.data.filename, header = TRUE, sep="\t")

    timepts.names <- input.data[1:num.timepts, 1]

    ## Remove first col i.e. the time point names
    input.data <- input.data[, -1]

  } else if (input.data.filename.ext[base::length(input.data.filename.ext)] == 'RData') {
    ## Loads an object named input.data
    base::load(input.data.filename)

    timepts.names <- 1:num.timepts
  }

  ## Begin: Replace original node names with {v1, v2, ..., vV}
  ## V = total number of nodes in the input data.
  ## The replacement is crucial as unexpected characters in the node names
  ## may generate an error in further computation.

  ## Save the original node names in case they are required in future
  orig.node.names <- base::colnames(input.data)

  node.names <- base::c()
  for (col.idx in 1:base::ncol(input.data))
  {
    new.node.name <- base::paste('v', base::as.character(col.idx), sep = '')
    node.names <- base::c(node.names, new.node.name)
  }
  base::rm(col.idx)
  base::colnames(input.data) <- node.names
  ## End: Replace original node names with {v1, v2, ..., vV}
  ## V = total number of node in the input data.

  num.nodes <- base::ncol(input.data)

  ## Max fanin must be restricted to 14. Because, it is empirically observed that bnstruct::learn.network()
  ## function can learn a BN with upto 15 nodes without segmentation fault, given a 32 GB main memory. A max
  ## fanin restriction of 14 limits the number of nodes in the to-be-learnt BN to 15 (1 target node and a
  ## maximum of 14 candidate parents).
  max.fanin <- base::min(num.nodes, 14)

  num.samples.per.timept <- (base::nrow(input.data) / num.timepts)

  ## If input data is already discretized
  if (is.discrete)
  {
    input.data.discr <- input.data

  } else
  {
    if (discr.algo == '')
    {
      stop('Please specify the value of discr.algo.')

    } else if (discr.algo == 'discretizeData.2L.wt.l')
    {
      input.data.discr <- discretizeData.2L.wt.l(input.data, input.wt.data.filename)

    } else if (discr.algo == 'discretizeData.2L.Tesla')
    {
      input.data.discr <- discretizeData.2L.Tesla(input.data)
    }

    base::save(input.data.discr, file = base::paste(output.dirname, 'input.data.discr.RData', sep = '/'))
  }
  ##------------------------------------------------------------
  ## End: Read input data
  ##------------------------------------------------------------

  ##------------------------------------------------------------
  ## Begin: Create a 3D array using discretized input data.
  ## Here,
  ## dim1 = time points,
  ## dim2 = nodes,
  ## dim3 = time series.
  ## It is useful for some CLR algos and the BN struct learning
  ## algos.
  ##------------------------------------------------------------
  input.data.discr.matrix <- base::data.matrix(input.data.discr)

  input.data.discr.3D <- base::array(NA, base::c(num.timepts, num.nodes, num.samples.per.timept),
                                     dimnames = base::c(base::list(timepts.names), base::list(node.names),
                                                        base::list(1:num.samples.per.timept)))

  for (sample.idx in 1:num.samples.per.timept) {
    start.row.idx <- (1 + (num.timepts * (sample.idx - 1)))
    end.row.idx <- (num.timepts * sample.idx)
    input.data.discr.3D[ , , sample.idx] <- input.data.discr.matrix[start.row.idx:end.row.idx, ]
  }
  base::rm(sample.idx)

  base::rm(input.data.discr.matrix)
  ##------------------------------------------------------------
  ## End: Create a 3D array using discretized input data.
  ##------------------------------------------------------------

  ##------------------------------------------------------------
  ## Begin: Learn MI Net Structure
  ##------------------------------------------------------------

  start.time <- base::proc.time() # start the timer

  ##------------------------------------------------------------
  ## Begin: Initialize mutual information network
  ##------------------------------------------------------------

  ## Initialize mutual info matrix for static CLR net
  mi.net.adj.matrix <- NULL

  ## Initialize mutual info matrix for time-varying CLR nets
  mi.net.adj.matrix.list <- NULL

  mut.info.matrix <- NULL
  if (clr.algo == 'CLR') {

    # Initialize mutual information matrix with zeroes
    mut.info.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                              dimnames = base::c(base::list(node.names), base::list(node.names)))

    if (mi.estimator == 'mi.pca.cmi') {
      ## Build mutual information matrix
      for (col.idx in 1:(num.nodes - 1)) {
        for (col.idx.2 in (col.idx + 1):num.nodes) {

          ## 'compute_cmi.R'
          mut.info <- ComputeCmiPcaCmi(input.data.discr[, col.idx], input.data.discr[, col.idx.2])

          mut.info.matrix[col.idx, col.idx.2] <- mut.info
          mut.info.matrix[col.idx.2, col.idx] <- mut.info
        }
        base::rm(col.idx.2)
      }
      base::rm(col.idx)

    } else if (mi.estimator == 'mi.empirical') {
      mut.info.matrix <- minet::build.mim(input.data.discr,
                                          estimator = 'mi.empirical',
                                          disc = 'none')

    } else if (mi.estimator == 'mi.mm') {
      mut.info.matrix <- minet::build.mim(input.data.discr,
                                          estimator = 'mi.mm',
                                          disc = 'none')
    }

    if (apply.aracne == TRUE) {

      mut.info.matrix.pre.aracne <- mut.info.matrix

      mut.info.matrix <- minet::aracne(mut.info.matrix)

      mut.info.matrix.post.aracne <- mut.info.matrix

      elapsed.time <- (base::proc.time() - start.time)
      base::writeLines('elapsed.time just after the ARACNE step= \n')
      base::print(elapsed.time)
      base::rm(elapsed.time)

      base::save(mut.info.matrix.pre.aracne,
           file = base::paste(output.dirname, 'mut.info.matrix.pre.aracne.RData', sep = '/'))

      base::save(mut.info.matrix.post.aracne,
           file = base::paste(output.dirname, 'mut.info.matrix.post.aracne.RData', sep = '/'))

      base::rm(mut.info.matrix.pre.aracne, mut.info.matrix.post.aracne)

    } else {
      ## apply.aracne == FALSE

      # writeLines('mut.info.matrix = \n')
      # print(mut.info.matrix)
      base::save(mut.info.matrix, file = base::paste(output.dirname, 'mut.info.matrix.RData', sep = '/'))
    }
  } else {
    ## clr.algo != 'CLR'

    base::rm(mut.info.matrix)
  }

  if ((clr.algo == 'CLR') | (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1') | (clr.algo == 'spearman')) {

    ## CLR net is not time-varying
    base::rm(mi.net.adj.matrix.list)

    mi.net.adj.matrix <- base::matrix(0, nrow = num.nodes, ncol = num.nodes,
                                dimnames = base::c(base::list(node.names), base::list(node.names)))

  } else if (clr.algo == 'CLR3') {

    ## CLR net is not static
    base::rm(mi.net.adj.matrix)

    ## Pre-allocate an empty list of length = number of time intervals
    num.time.ivals <- (num.timepts - 1)
    mi.net.adj.matrix.list <- base::vector(mode = 'list', length = num.time.ivals)
    base::rm(num.time.ivals)
  }
  ##------------------------------------------------------------
  ## End: Initialize mutual information network
  ##------------------------------------------------------------

  # entropy.matrix <- computEntropy(input.data.discr) #----Verify the name

  ## Initialize filename where 'mi.net.adj.matrix.list' is to be saved
  ## in case 'clr.algo == CLR3'
  mi.net.adj.matrix.list.filename <- NULL
  if ((clr.algo == 'CLR') | (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1')) {
    base::rm(mi.net.adj.matrix.list.filename)
  }

  ## source('learn_mi_net_struct.R')
  if (clr.algo == 'CLR') {
    # mi.net.adj.matrix <- LearnMiNetStructZstat(mut.info.matrix, mi.net.adj.matrix, entropy.matrix, alpha)
    # mi.net.adj.matrix <- LearnMiNetStructClr(mut.info.matrix, mi.net.adj.matrix, num.nodes)
    mi.net.adj.matrix <- LearnClrNetMfi(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin, output.dirname)

  } else if (clr.algo == 'CLR2') {
    mi.net.adj.matrix <- LearnClr2NetMfi(input.data.discr, num.nodes, node.names, num.timepts,
                                         max.fanin, output.dirname, mi.net.adj.matrix)

  } else if (clr.algo == 'CLR2.1') {
    mi.net.adj.matrix <- LearnClrNetMfiVer2.1(input.data.discr, num.nodes, node.names, num.timepts,
                                              max.fanin, output.dirname, mi.net.adj.matrix)

  } else if (clr.algo == 'CLR3') {
    mi.net.adj.matrix.list <- LearnClr3NetMfi(input.data.discr.3D, num.nodes, node.names, num.timepts,
                                              max.fanin, mi.net.adj.matrix.list)

    ## Since 'mi.net.adj.matrix.list' is very large, save it in a specific file
    ## and remove it. Then load it when necessary. No need to retain it in the
    ## workspace even when not required.
    mi.net.adj.matrix.list.filename <- base::paste(output.dirname, 'mi.net.adj.matrix.list.RData', sep = '/')
    base::save(mi.net.adj.matrix.list, file = mi.net.adj.matrix.list.filename)
    base::rm(mi.net.adj.matrix.list)
  }

  elapsed.time <- (base::proc.time() - start.time) # Check time taken by CLR
  base::writeLines('elapsed.time just after the CLR step= \n')
  base::print(elapsed.time)
  base::rm(elapsed.time)

  if (clr.algo == 'CLR') {
    base::rm(mut.info.matrix)
  }

  if ((clr.algo == 'CLR') | (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1')) {
    # writeLines('\n mi.net.adj.matrix = \n')
    # print(mi.net.adj.matrix)
    base::save(mi.net.adj.matrix, file = base::paste(output.dirname, 'mi.net.adj.matrix.RData', sep = '/'))
    ## Check max number of nbrs for a node in the MI net
    # writeLines('\n Max number of nbrs = \n')
    # print(max(colSums(mi.net.adj.matrix)))

    ## Identify which nodes have the max number of nbrs
    # writeLines('\n Nodes having max number of nbrs = \n')
    # print(colnames(mi.net.adj.matrix[, colSums(mi.net.adj.matrix) == max(colSums(mi.net.adj.matrix))]))
  }

  # stop('MI net struct learning completed.')

  # Rgraphviz::plot(as(mi.net.adj.matrix, 'graphNEL'))
  ##------------------------------------------------------------
  ## End: Learn MI Net Structure
  ##------------------------------------------------------------

  ##------------------------------------------------------------
  ## Begin: Learn Network Structures
  ##------------------------------------------------------------
  unrolled.DBN.adj.matrix.list <- NULL

  if ((clr.algo == 'CLR') | (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1')) {
    unrolled.DBN.adj.matrix.list <- learnDbnStructMo1Layer3dParDeg1_v2(input.data.discr.3D, mi.net.adj.matrix,
                                                                       num.discr.levels, num.nodes, num.timepts,
                                                                       max.fanin, node.names, clr.algo)
    base::rm(mi.net.adj.matrix)

  } else if (clr.algo == 'CLR3') {

    num.time.ivals <- (num.timepts - 1)
    unrolled.DBN.adj.matrix.list <- base::vector(mode = 'list', length = num.time.ivals)
    # rm(num.time.ivals)

    ## Adjacency matrix for a time-interval-specific DBN
    time.ival.spec.dbn.adj.matrix <- base::matrix(0, nrow = num.nodes,
                                            ncol = num.nodes,
                                            dimnames = base::c(base::list(node.names),
                                                              base::list(node.names)))
    for (time.ival.idx in 1:num.time.ivals) {
      unrolled.DBN.adj.matrix.list[[time.ival.idx]] <- time.ival.spec.dbn.adj.matrix
    }
    base::rm(time.ival.idx)

    base::rm(num.time.ivals, time.ival.spec.dbn.adj.matrix)

    unrolled.DBN.adj.matrix.list <- LearnDbnStructMo1Clr3Ser(input.data.discr.3D, mi.net.adj.matrix.list.filename,
                                                             num.discr.levels, num.nodes, num.timepts, max.fanin,
                                                             node.names, unrolled.DBN.adj.matrix.list)
    base::rm(mi.net.adj.matrix.list.filename)
  }

  base::save(unrolled.DBN.adj.matrix.list, file = base::paste(output.dirname, 'unrolled.DBN.adj.matrix.list.RData', sep = '/'))
  base::rm(input.data.discr.3D)

  ## Learn the rolled DBN adj matrix
  ## source(paste(init.path, 'rollDbn.R', sep = '/'))
  # rolled.DBN.adj.matrix <- rollDbn(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, roll.method, allow.self.loop)
  # rolled.DBN.adj.matrix <- rollDbn(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, 'any', FALSE)
  rolled.DBN.adj.matrix <- rollDbn_v2(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix.list,
                                      'any', allow.self.loop)
  di.net.adj.matrix <- rolled.DBN.adj.matrix
  base::rm(rolled.DBN.adj.matrix)

  # writeLines('\n di.net.adj.matrix = \n')
  # print(di.net.adj.matrix)
  ## Change the node names back to the original node names
  base::rownames(di.net.adj.matrix) <- orig.node.names
  base::colnames(di.net.adj.matrix) <- orig.node.names
  base::save(di.net.adj.matrix, file = base::paste(output.dirname, 'di.net.adj.matrix.RData', sep = '/'))

  ## Create an '.sif' file equivalent to the directed net adjacency matrix
  ## that is readable in Cytoscape.
  adjmxToSif(di.net.adj.matrix, output.dirname)
  # rm(unrolled.DBN.adj.matrix)
  base::rm(unrolled.DBN.adj.matrix.list)
  ##------------------------------------------------------------
  ## End: Learn Network Structures
  ##------------------------------------------------------------

  ## If the true rolled network is known a prior, then evaluate the performance
  ## metrics of the predicted rolled network.
  if (true.net.filename != '')
  {
    predicted.net.adj.matrix <- di.net.adj.matrix

    ## Loads R obj 'true.net.adj.matrix'
    true.net.adj.matrix <- NULL
    base::load(true.net.filename)


    ## Begin: Create the format for result
    Result <- base::matrix(0, nrow = 1, ncol = 11)
    base::colnames(Result) <- base::list('TP', 'TN', 'FP', 'FN', 'TPR', 'FPR', 'FDR', 'PPV', 'ACC', 'MCC',  'F')
    # ## End: Create the format for result

    ResultVsTrue <- calcPerfDiNet(predicted.net.adj.matrix, true.net.adj.matrix, Result, num.nodes)
    base::rm(Result)
    base::writeLines('Result TGS vs True = \n')
    base::print(ResultVsTrue)
    base::rm(ResultVsTrue)
  }

  base::rm(di.net.adj.matrix)

  elapsed.time <- (base::proc.time() - start.time) # Stop the timer
  base::writeLines('elapsed.time = \n')
  base::print(elapsed.time)
  base::rm(elapsed.time)

  base::sink()
  base::close(output.file.conn)
  ##------------------------------------------------------------
  ## Begin: Save R session info in a File
  ##-----------------------------------------------------------
  session.file <- base::file(base::paste(output.dirname, 'sessionInfo.txt', sep = '/'), open = "wt")
  base::sink(session.file)
  utils::sessionInfo()
  base::sink()
  base::close(session.file)
  ##------------------------------------------------------------
  ## End: Save R session info in a File
  ##------------------------------------------------------------

  ##------------------------------------------------------------
  ## End: Main program
  ##------------------------------------------------------------

}
