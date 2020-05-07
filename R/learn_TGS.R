#'Implement the TGS Algorithm
#'
#'The TGS algorithm takes a time-series gene expression dataset as input. It
#'analyses the data and reconstructs the underlying temporal sequence of gene
#'regulatory events. The reconstructed output is given in the form of
#'time-varying gene regulatory networks (GRNs). The TGS algorithm is extremely
#'time-efficient and hence suitable for processing large datasets with hundreds
#'to thousands of genes. More details about the algorithm can be found at
#'Saptarshi Pyne, Alok Ranjan Kumar, and Ashish Anand. Rapid reconstruction of
#'time-varying gene regulatory networks. IEEE/ACM Transactions on Computational
#'Biology and Bioinformatics, 17(1):278--291, Jan--Feb 2020.
#'
#'@param isfile Numeric. 1 or 0. 1 if input arguments are given in a json file.
#'  Otherwise, 0.
#'
#'@param json.file Character string. Absolute path to the JSON file if
#'  \code{isfile = 1}.
#'
#'@param input.dirname Character string. Absolute path to the directory where
#'  input files are kept. By default, the current working directory.
#'
#'@param input.data.filename Character string. Name of the file containing the
#'  input data. It can either be a '.tsv' file or an '.RData' file. \itemize{
#'  \item If it is a '.tsv' file, the first column should have the time point
#'  IDs. The only exception is the (1, 1)-th cell. This cell should be reserved
#'  for the column header. The column header can be anything, such as - 'Time'.
#'  The first row, excluding the (1, 1)-th cell, must have the gene names. The
#'  rest of the cells should contain the corresponding value. For example, the
#'  cell with row ID 't1' and column ID 'G2' would represent the expression
#'  value of gene 'G2' at time point 't1'.
#'
#'  \item If it is an '.RData' file, the underlying object should be a matrix
#'  named 'input.data'. In 'input.data', the row names must represent the time
#'  point IDs. On the other hand, the column names must represent the gene
#'  names. Therefore, the (i, j)-th cell of the matrix contains the expression
#'  value of the j-th gene at the i-th time point. }
#'
#'  For both '.tsv' and '.RData' input, multiple rows with the same time point
#'  ID represent multiple replicates at the same time point. In other words,
#'  these rows belong to the same time point in different time series. The time
#'  points belonging to the same time series must be together and in ascending
#'  order. An exemplary dataset with three genes \{G1, G2, G3\}, two time points
#'  \{t1, t2\} and two time series is shown below.
#'
#'  \tabular{rrrr}{ \strong{Time} \tab \strong{G1} \tab \strong{G2} \tab
#'  \strong{G3} \cr t1 \tab 0.8272480342 \tab 0.7257430901 \tab 0.3894130418 \cr
#'  t2 \tab 0.6542518342 \tab 0.6470658823 \tab 0.5088904888 \cr t1 \tab
#'  0.3519554463 \tab 0.3551279726 \tab 0.3207993604 \cr t2 \tab 0.4871730974
#'  \tab 0.3706990326 \tab 0.447523615 }
#'
#'@param num.timepts Numeric. Positive integer greater than 1. Number of
#'  distinct time points.
#'
#'@param true.net.filename Character string. Name of the file containing the
#'  true network. In case it is non-empty, the name should refer to an '.RData'
#'  file. The '.RData' file must have an object named 'true.net.adj.matrix'. The
#'  object can either be a matrix or a list. \itemize{ \item If the object is a
#'  matrix, then it represents the true summary GRN. The row names and column
#'  names should be the gene names. Each cell can contain a value of 1 or 0. If
#'  the (i, j)-th cell contains 1, then there exists an edge from the i-th gene
#'  to the j-th gene. Otherwise, the edge does not exist. An example with three
#'  genes \{G1, G2, G3\} is given below. \tabular{rrrr}{ \tab \strong{G1} \tab
#'  \strong{G2} \tab \strong{G3} \cr \strong{G1} \tab 0 \tab 1 \tab 0 \cr
#'  \strong{G2} \tab 0 \tab 0 \tab 0 \cr \strong{G3} \tab 1 \tab 0 \tab 0 }
#'
#'  \item If the object is a list, then it represents the true time-varying
#'  GRNs. The length of the list must be equal to the number of time intervals,
#'  which is \code{(num.timepts - 1)}. Each element in the list should be a
#'  matrix. The p-th matrix represents the true GRN corresponding to the p-th
#'  time interval. The row names and column names should be the gene names. Each
#'  cell can contain a value of 1 or 0. If the (i, j)-th cell contains 1, then
#'  there exists an edge from the i-th gene to the j-th gene. Otherwise, the
#'  edge does not exist. }
#'
#'@param input.wt.data.filename Character string. Name of the file containing
#'  the Wild Type expressions of the genes. If non-empty, then must be a '.tsv'
#'  file. The first row should contain the names of the genes. Only exception is
#'  the (1, 1)-th cell which should be empty. The second row should have the
#'  wild type expressions. Therefore, the (2, j)-th cell must contain the wild
#'  type expression of the j-th gene. Again the only exception is the (2, 1)-th
#'  cell which should be empty. An example with three genes \{G1, G2, G3\} is
#'  given below.
#'
#'  \tabular{rrrr}{ <empty> \tab \strong{G1} \tab \strong{G2} \tab \strong{G3}
#'  \cr <empty>  \tab 0.5298713 \tab 0.5174261 \tab 0.8181522 }
#'
#'
#'@param is.discrete Logical. TRUE or FALSE. TRUE if the input data is discrete.
#'  Otherwise, FALSE.
#'
#'@param num.discr.levels Numeric. Positive integer greater than 1. Number of
#'  discrete levels that each gene has (if the input data is discrete) or each
#'  gene should have (if the input data needs to be discretised).
#'
#'@param discr.algo Character string. Name of the discretisation algorithm to be
#'  used when the input data needs to be discretised. The available algorithms
#'  are -- 'discretizeData.2L.Tesla' and 'discretizeData.2L.wt.l'. If you
#'  choose algorithm 'discretizeData.2L.wt.l', please provide the wild type
#'  data using argument \code{input.wt.data.filename}.
#'
#'@param mi.estimator Character string. Name of the algorithm for estimating
#'  mutual informations. There is only one algorithm available at this moment.
#'  It is 'mi.pca.cmi'.
#'
#'@param apply.aracne Logical. TRUE or FALSE. TRUE if you wish to apply ARACNE
#'  for refining the mutual information matrix. Otherwise, FALSE.
#'
#'@param clr.algo Character string. Name of the context likelihood relatedness
#'  (CLR) algorithm to use. The available algorithms are -- 'CLR', 'CLR2',
#'  'CLR2.1', 'CLR3' and 'spearman'.
#'
#'@param max.fanin Numeric. Positive integer. Maximum number of regulators each
#'  gene can have.
#'
#'@param allow.self.loop Logical. TRUE or FALSE. TRUE if you wish to allow self
#'  loops. Otherwise, FALSE.
#'
#'@param scoring.func Character string. Name of the scoring function to use. At
#'  this moment, the only available option is 'BIC'.
#'
#'@param output.dirname Character string. File path to a directory where output
#'  files are to be saved. There are three options. \emph{Option 1:} It can be
#'  the absolute path to an existing directory. \emph{Option 2:} It can also be
#'  the absolute path to a non-existing directory. In this case, the directory
#'  will be created. \emph{Option 3 (default):} If provided an empty string,
#'  then it will be the current working directory.
#'
#'@details The function does not return any values. Instead, it outputs a set of
#'  files and saves them under the directory specified by \code{output.dirname}.
#'  The output files are described in Section 'Value'.
#'
#'@return \describe{ \item{input.data.discr.RData}{ Discretised version of the
#'  input data. This file is created only if the input data is not discretised
#'  as specified by input argument 'is.discrete'. }
#'
#'  \item{mut.info.matrix.RData}{ Mutual information matrix of the given genes.
#'  This RData file contains a matrix named 'mut.info.matrix'. The (i, j)-th
#'  cell of the matrix represents the mutual information between the i-th and
#'  j-th genes. This is a symmetric matrix. }
#'
#'  \item{mi.net.adj.matrix.wt.RData}{ Weighted Mutual information network of
#'  the given genes. This RData file contains a matrix named
#'  'mi.net.adj.matrix.wt'. The (i, j)-th cell of the matrix represents the
#'  weight of the edge from the i-th gene to the j-th gene. The edge weight is a
#'  non-negative real number. }
#'
#'  \item{mi.net.adj.matrix.RData}{ Unweighted Mutual information network of the
#'  given genes. This RData file contains a matrix named 'mi.net.adj.matrix'.
#'  Each cell of the matrix contains a value of 1 or 0. If the (i, j)-th cell
#'  contains 1, then there exists an edge from the i-th gene to the j-th gene.
#'  Otherwise, the edge does not exist. }
#'
#'  \item{unrolled.DBN.adj.matrix.list.RData}{ Reconstructed time-varying GRNs.
#'  This RData file contains a list named 'unrolled.DBN.adj.matrix.list'. The
#'  length of the list is equal to the total number of time intervals, which is
#'  \code{(num.timepts - 1)}. Each element in the list is a network adjacency
#'  matrix. The p-th element in the list represents the adjacency matrix of the
#'  GRN corresponding to the p-th time interval. In this adjacency matrix, each
#'  cell contains a value of 1 or 0. If the (i, j)-th cell contains 1, then
#'  there exists a directed edge from the i-th gene to the j-th gene. Otherwise,
#'  the edge does not exist. }
#'
#'  \item{di.net.adj.matrix.RData}{ Rolled GRN. This RData file contains a
#'  matrix named 'di.net.adj.matrix'. Each cell in the matrix contains a value
#'  of 1 or 0. If the (i, j)-th cell contains 1, then there exists an edge from
#'  the i-th gene to the j-th gene. Otherwise, the edge does not exist. }
#'
#'  \item{net.sif}{ Rolled GRN in the SIF format compatible with Cytoscape. }
#'
#'  \item{Result.RData}{ Correctness metrics. This file is created only if true
#'  network is given through input argument 'true.net.filename'. Inside this
#'  RData file, there is a matrix named 'Result'. The columns represent the
#'  correctness metrics, such as - TP (number of true positive predictions) and
#'  FP (number of false positive predictions). The rows depend upon the nature
#'  of the true network. If the true network is time-varying GRNs, then the
#'  number of rows is equal to the number of time intervals. In that case, the
#'  p-th row contains the correctness metrics of the reconstructed GRN
#'  corresponding to the p-th time interval. On the other hand, if the true
#'  network is a summary GRN, then there exists only one row. This row
#'  represents the correctness metrics of the rolled GRN. }
#'
#'  \item{output.txt}{ Console output. }
#'
#'  \item{sessionInfo.txt}{ R session information. } }
#'
#' @examples
#' \dontrun{
#'   TGS::LearnTgs(
#'   isfile = 0,
#'   json.file = '',
#'   input.dirname = 'C:/GitHub/TGS/inst/extdata',
#'   input.data.filename = 'InSilicoSize10-Yeast1-trajectories.tsv',
#'   num.timepts = 21,
#'   true.net.filename = 'DREAM3GoldStandard_InSilicoSize10_Yeast1_TrueNet.RData',
#'   input.wt.data.filename = 'InSilicoSize10-Yeast1-null-mutants.tsv',
#'   is.discrete = FALSE,
#'   num.discr.levels = 2,
#'   discr.algo = 'discretizeData.2L.wt.l',
#'   mi.estimator = 'mi.pca.cmi',
#'   apply.aracne = FALSE,
#'   clr.algo = 'CLR',
#'   max.fanin = 14,
#'   allow.self.loop = FALSE,
#'   scoring.func = 'BIC',
#'   output.dirname = 'C:/GitHub/TGS/inst/extdata/Output_Ds10n')
#'
#'   TGS::LearnTgs(
#'   isfile = 0,
#'   json.file = '',
#'   input.dirname = 'C:/GitHub/TGS/inst/extdata',
#'   input.data.filename = 'edi-data-10n.tsv',
#'   num.timepts = 21,
#'   true.net.filename = 'edi.net.10.adj.mx.RData',
#'   input.wt.data.filename = '',
#'   is.discrete = FALSE,
#'   num.discr.levels = 2,
#'   discr.algo = 'discretizeData.2L.Tesla',
#'   mi.estimator = 'mi.pca.cmi',
#'   apply.aracne = FALSE,
#'   clr.algo = 'CLR',
#'   max.fanin = 14,
#'   allow.self.loop = TRUE,
#'   scoring.func = 'BIC',
#'   output.dirname = 'C:/GitHub/TGS/inst/extdata/Output_Ed10n')
#' }
#'
#'@export
LearnTgs <- function(isfile = 0,
                     json.file = '',
                     input.dirname = '',
                     input.data.filename = '',
                     num.timepts = 2,
                     true.net.filename = '',
                     input.wt.data.filename = '',
                     is.discrete = TRUE,
                     num.discr.levels = 2,
                     discr.algo = '',
                     mi.estimator = 'mi.pca.cmi',
                     apply.aracne = FALSE,
                     clr.algo = 'CLR',
                     max.fanin = 14,
                     allow.self.loop = TRUE,
                     scoring.func = 'BIC',
                     output.dirname = '') {
  if (isfile != 0) {
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
    scoring.func <- input.params$scoring.func
    input.dirname <- input.params$input.dirname
    output.dirname <- input.params$output.dirname
    rm(input.params)

    ##------------------------------------------------------------
    ## End: Read User-defined input Params
    ##------------------------------------------------------------
  }

  if (input.dirname == '') {
    input.dirname <- base::getwd()
  }

  if (output.dirname == '') {
    output.dirname <- base::getwd()

  } else if (!base::file.exists(output.dirname)) {
    ##------------------------------------------------------------
    ## Begin: Create the output directory
    ##------------------------------------------------------------
    if (base::.Platform$OS.type == 'windows') {
      ## Convert directory path to canonical form for Windows OS.
      ## It raises the warning if the directory does not exist, which
      ## is expected. Therefore, please ignore the warning.
      output.dirname <-
        base::normalizePath(output.dirname,
                            winslash = '\\',
                            mustWork = NA)

      shell(
        base::paste('mkdir ', output.dirname, sep = ''),
        intern = TRUE,
        mustWork = TRUE
      )

    } else if (base::.Platform$OS.type == 'unix') {
      base::system(base::paste('mkdir ', output.dirname, sep = ''))
    }
    ##------------------------------------------------------------
    ## End: Create the output directory
    ##------------------------------------------------------------
  }

  input.data.filename <-
    base::paste(input.dirname, input.data.filename, sep = '/')

  if (true.net.filename != '') {
    true.net.filename <-
      base::paste(input.dirname, true.net.filename, sep = '/')
  }

  if (input.wt.data.filename != '') {
    input.wt.data.filename <-
      base::paste(input.dirname, input.wt.data.filename, sep = '/')
  }

  ##------------------------------------------------------------
  ## Begin: Main program
  ##------------------------------------------------------------

  base::print('The output directory name is:')
  base::print(output.dirname)
  base::print('') ## to append a blank line

  ## Save console output in a file named 'output.txt' inside the output directory.
  output.filename <-
    base::paste(output.dirname, 'output.txt', sep = '/')
  output.file.conn <- base::file(output.filename, open = "wt")
  base::sink(output.file.conn)

  ##------------------------------------------------------------
  ## Begin: Read input data file
  ##------------------------------------------------------------

  ## Begin: Find file extension of the input data file. Only '.tsv' and '.RData'
  ## are allowed.
  ## Split the string at every '.' and consider the last substring as the
  ## file extension.
  input.data.filename.ext <-
    base::unlist(base::strsplit(input.data.filename, '[.]'))
  ## End: Find file extension of the input data file. Only '.tsv' and '.RData'
  ## are allowed.

  ## Initialize input data
  input.data <- NULL
  if (input.data.filename.ext[base::length(input.data.filename.ext)] == 'tsv') {
    input.data <-
      utils::read.table(input.data.filename, header = TRUE, sep = "\t")

    timepts.names <- input.data[1:num.timepts, 1]

    ## Remove first col i.e. the time point names
    input.data <- input.data[,-1]

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
    new.node.name <-
      base::paste('v', base::as.character(col.idx), sep = '')
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

  } else {
    if (discr.algo == '') {
      stop('Please specify the value of discr.algo.')

    } else if (discr.algo == 'discretizeData.2L.wt.l') {
      input.data.discr <-
        discretizeData.2L.wt.l(input.data, input.wt.data.filename)

    } else if (discr.algo == 'discretizeData.2L.Tesla') {
      input.data.discr <-
        discretizeData.2L.Tesla(input.data)
    }

    base::save(
      input.data.discr,
      file = base::paste(output.dirname, 'input.data.discr.RData', sep = '/')
    )
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

  input.data.discr.3D <-
    base::array(
      NA,
      base::c(num.timepts, num.nodes, num.samples.per.timept),
      dimnames = base::c(
        base::list(timepts.names),
        base::list(node.names),
        base::list(1:num.samples.per.timept)
      )
    )

  for (sample.idx in 1:num.samples.per.timept) {
    start.row.idx <- (1 + (num.timepts * (sample.idx - 1)))
    end.row.idx <- (num.timepts * sample.idx)
    input.data.discr.3D[, , sample.idx] <-
      input.data.discr.matrix[start.row.idx:end.row.idx,]
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
    mut.info.matrix <-
      base::matrix(
        0,
        nrow = num.nodes,
        ncol = num.nodes,
        dimnames = base::c(base::list(node.names), base::list(node.names))
      )

    if (mi.estimator == 'mi.pca.cmi') {
      ## Build mutual information matrix
      for (col.idx in 1:(num.nodes - 1)) {
        for (col.idx.2 in (col.idx + 1):num.nodes) {
          ## 'compute_cmi.R'
          mut.info <-
            ComputeCmiPcaCmi(input.data.discr[, col.idx],
                                  input.data.discr[, col.idx.2])

          mut.info.matrix[col.idx, col.idx.2] <- mut.info
          mut.info.matrix[col.idx.2, col.idx] <- mut.info
        }
        base::rm(col.idx.2)
      }
      base::rm(col.idx)

    } else if (mi.estimator == 'mi.empirical') {
      mut.info.matrix <-
        minet::build.mim(input.data.discr,
                         estimator = 'mi.empirical',
                         disc = 'none')

    } else if (mi.estimator == 'mi.mm') {
      mut.info.matrix <-
        minet::build.mim(input.data.discr,
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

      base::save(
        mut.info.matrix.pre.aracne,
        file = base::paste(
          output.dirname,
          'mut.info.matrix.pre.aracne.RData',
          sep = '/'
        )
      )

      base::save(
        mut.info.matrix.post.aracne,
        file = base::paste(
          output.dirname,
          'mut.info.matrix.post.aracne.RData',
          sep = '/'
        )
      )

      base::rm(mut.info.matrix.pre.aracne,
               mut.info.matrix.post.aracne)

    } else {
      ## apply.aracne == FALSE

      # writeLines('mut.info.matrix = \n')
      # print(mut.info.matrix)
      base::save(
        mut.info.matrix,
        file = base::paste(output.dirname, 'mut.info.matrix.RData', sep = '/')
      )
    }
  } else {
    ## clr.algo != 'CLR'

    base::rm(mut.info.matrix)
  }

  if ((clr.algo == 'CLR') |
      (clr.algo == 'CLR2') |
      (clr.algo == 'CLR2.1') | (clr.algo == 'spearman')) {
    ## CLR net is not time-varying
    base::rm(mi.net.adj.matrix.list)

    mi.net.adj.matrix <-
      base::matrix(
        0,
        nrow = num.nodes,
        ncol = num.nodes,
        dimnames = base::c(base::list(node.names), base::list(node.names))
      )

  } else if (clr.algo == 'CLR3') {
    ## CLR net is not static
    base::rm(mi.net.adj.matrix)

    ## Pre-allocate an empty list of length = number of time intervals
    num.time.ivals <- (num.timepts - 1)
    mi.net.adj.matrix.list <-
      base::vector(mode = 'list', length = num.time.ivals)
    base::rm(num.time.ivals)
  }
  ##------------------------------------------------------------
  ## End: Initialize mutual information network
  ##------------------------------------------------------------

  # entropy.matrix <- computEntropy(input.data.discr) #----Verify the name

  ## Initialize filename where 'mi.net.adj.matrix.list' is to be saved
  ## in case 'clr.algo == CLR3'
  mi.net.adj.matrix.list.filename <- NULL
  if ((clr.algo == 'CLR') |
      (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1')) {
    base::rm(mi.net.adj.matrix.list.filename)
  }

  ## source('learn_mi_net_struct.R')
  if (clr.algo == 'CLR') {
    # mi.net.adj.matrix <- LearnMiNetStructZstat(mut.info.matrix, mi.net.adj.matrix, entropy.matrix, alpha)
    # mi.net.adj.matrix <- LearnMiNetStructClr(mut.info.matrix, mi.net.adj.matrix, num.nodes)
    mi.net.adj.matrix <-
      LearnClrNetMfi(mut.info.matrix,
                          mi.net.adj.matrix,
                          num.nodes,
                          max.fanin,
                          output.dirname)

  } else if (clr.algo == 'CLR2') {
    mi.net.adj.matrix <-
      LearnClr2NetMfi(
        input.data.discr,
        num.nodes,
        node.names,
        num.timepts,
        max.fanin,
        output.dirname,
        mi.net.adj.matrix
      )

  } else if (clr.algo == 'CLR2.1') {
    mi.net.adj.matrix <-
      LearnClrNetMfiVer2.1(
        input.data.discr,
        num.nodes,
        node.names,
        num.timepts,
        max.fanin,
        output.dirname,
        mi.net.adj.matrix
      )

  } else if (clr.algo == 'CLR3') {
    mi.net.adj.matrix.list <-
      LearnClr3NetMfi(
        input.data.discr.3D,
        num.nodes,
        node.names,
        num.timepts,
        max.fanin,
        mi.net.adj.matrix.list
      )

    ## Since 'mi.net.adj.matrix.list' is very large, save it in a specific file
    ## and remove it. Then load it when necessary. No need to retain it in the
    ## workspace even when not required.
    mi.net.adj.matrix.list.filename <-
      base::paste(output.dirname, 'mi.net.adj.matrix.list.RData', sep = '/')
    base::save(mi.net.adj.matrix.list, file = mi.net.adj.matrix.list.filename)
    base::rm(mi.net.adj.matrix.list)
  }

  elapsed.time <-
    (base::proc.time() - start.time) # Check time taken by CLR
  base::writeLines('elapsed.time just after the CLR step= \n')
  base::print(elapsed.time)
  base::rm(elapsed.time)

  if (clr.algo == 'CLR') {
    base::rm(mut.info.matrix)
  }

  if ((clr.algo == 'CLR') |
      (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1')) {
    # writeLines('\n mi.net.adj.matrix = \n')
    # print(mi.net.adj.matrix)
    base::save(
      mi.net.adj.matrix,
      file = base::paste(output.dirname, 'mi.net.adj.matrix.RData', sep = '/')
    )
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

  if ((clr.algo == 'CLR') |
      (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1')) {
    unrolled.DBN.adj.matrix.list <-
      learnDbnStructMo1Layer3dParDeg1_v2(
        input.data.discr.3D,
        mi.net.adj.matrix,
        num.discr.levels,
        num.nodes,
        num.timepts,
        max.fanin,
        node.names,
        clr.algo
      )

    base::rm(mi.net.adj.matrix)

  } else if (clr.algo == 'CLR3') {
    num.time.ivals <- (num.timepts - 1)
    unrolled.DBN.adj.matrix.list <-
      base::vector(mode = 'list', length = num.time.ivals)
    # rm(num.time.ivals)

    ## Adjacency matrix for a time-interval-specific DBN
    time.ival.spec.dbn.adj.matrix <-
      base::matrix(
        0,
        nrow = num.nodes,
        ncol = num.nodes,
        dimnames = base::c(base::list(node.names),
                           base::list(node.names))
      )
    for (time.ival.idx in 1:num.time.ivals) {
      unrolled.DBN.adj.matrix.list[[time.ival.idx]] <-
        time.ival.spec.dbn.adj.matrix
    }
    base::rm(time.ival.idx)

    base::rm(num.time.ivals, time.ival.spec.dbn.adj.matrix)

    unrolled.DBN.adj.matrix.list <-
      LearnDbnStructMo1Clr3Ser(
        input.data.discr.3D,
        mi.net.adj.matrix.list.filename,
        num.discr.levels,
        num.nodes,
        num.timepts,
        max.fanin,
        node.names,
        unrolled.DBN.adj.matrix.list
      )
    base::rm(mi.net.adj.matrix.list.filename)
  }

  base::save(
    unrolled.DBN.adj.matrix.list,
    file = base::paste(output.dirname, 'unrolled.DBN.adj.matrix.list.RData', sep = '/')
  )
  base::rm(input.data.discr.3D)

  ## Learn the rolled DBN adj matrix
  if (!base::is.null(unrolled.DBN.adj.matrix.list)) {
    ## 'unrolled.DBN.adj.matrix.list' is NULL if any only if
    ## none of the target nodes have any candidate parents in
    ## their corresponding shortlists

    # rolled.DBN.adj.matrix <- rollDbn(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, roll.method, allow.self.loop)
    # rolled.DBN.adj.matrix <- rollDbn(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, 'any', FALSE)
    rolled.DBN.adj.matrix <-
      rollDbn_v2(
        num.nodes,
        node.names,
        num.timepts,
        unrolled.DBN.adj.matrix.list,
        'any',
        allow.self.loop
      )
    di.net.adj.matrix <- rolled.DBN.adj.matrix
    base::rm(rolled.DBN.adj.matrix)

    # writeLines('\n di.net.adj.matrix = \n')
    # print(di.net.adj.matrix)
    ## Change the node names back to the original node names
    base::rownames(di.net.adj.matrix) <- orig.node.names
    base::colnames(di.net.adj.matrix) <- orig.node.names
    base::save(
      di.net.adj.matrix,
      file = base::paste(output.dirname, 'di.net.adj.matrix.RData', sep = '/')
    )

    ## Create an '.sif' file equivalent to the directed net adjacency matrix
    ## that is readable in Cytoscape.
    adjmxToSif(di.net.adj.matrix, output.dirname)
    # rm(unrolled.DBN.adj.matrix)
    # base::rm(unrolled.DBN.adj.matrix.list)
    ##------------------------------------------------------------
    ## End: Learn Network Structures
    ##------------------------------------------------------------


    ##------------------------------------------------------------
    ## Begin: Calc performance metrics if true net(s) is known
    ##------------------------------------------------------------
    if (true.net.filename != '') {
      ## Load R obj 'true.net.adj.matrix'
      true.net.adj.matrix <- NULL
      base::load(true.net.filename)

      ## Begin: Create the format for result
      Result <- base::matrix(0, nrow = 1, ncol = 11)
      base::colnames(Result) <-
        base::list('TP',
                   'TN',
                   'FP',
                   'FN',
                   'TPR',
                   'FPR',
                   'FDR',
                   'PPV',
                   'ACC',
                   'MCC',
                   'F1')
      ## End: Create the format for result

      if (base::is.matrix(true.net.adj.matrix)) {
        ## True net is time-invariant. Therefore,
        ## 'true.net.adj.matrix' is a single matrix.

        predicted.net.adj.matrix <- di.net.adj.matrix

        ResultVsTrue <-
          calcPerfDiNet(predicted.net.adj.matrix,
                             true.net.adj.matrix,
                             Result,
                             num.nodes)
        base::writeLines('Prediction vs Truth = \n')
        base::print(ResultVsTrue)
        base::rm(ResultVsTrue)

      } else if (base::is.list(true.net.adj.matrix)) {
        ## True nets are time-varying. Therefore,
        ## 'true.net.adj.matrix' is a list of matrices.

        for (net.idx in 1:base::length(unrolled.DBN.adj.matrix.list)) {
          predicted.net.adj.matrix <-
            unrolled.DBN.adj.matrix.list[[net.idx]]

          ResultVsTrue <-
            calcPerfDiNet(predicted.net.adj.matrix,
                               true.net.adj.matrix[[net.idx]],
                               Result,
                               num.nodes)
          Result <-
            base::rbind(Result,
                        base::matrix(
                          ResultVsTrue[1,],
                          nrow = 1,
                          ncol = ncol(Result)
                        ))

          # rm(ResultVsTrue)
        }
        base::rm(net.idx)

        ## Print mean performance averaged over
        ## all time-varying networks
        ResultVsTrue <- base::colMeans(Result)
        ResultVsTrue <-
          base::matrix(colMeans(Result),
                       nrow = 1,
                       ncol = ncol(Result))
        base::colnames(ResultVsTrue) <- base::colnames(Result)
        base::writeLines('Prediction vs Truth = \n')
        base::print(ResultVsTrue)
        base::rm(ResultVsTrue)
      }

      base::save(Result,
                 file =
                   base::paste(output.dirname, 'Result.RData', sep = '/'))
      base::rm(Result)
    }
    base::rm(di.net.adj.matrix)

  }
  base::rm(unrolled.DBN.adj.matrix.list)
  ##------------------------------------------------------------
  ## End: Calc performance metrics if true net(s) is known
  ##------------------------------------------------------------

  ## Stop the timer
  elapsed.time <- (base::proc.time() - start.time)
  base::writeLines('elapsed.time = \n')
  base::print(elapsed.time)
  base::rm(elapsed.time)

  base::sink()
  base::close(output.file.conn)
  ##------------------------------------------------------------
  ## Begin: Save R session info in a File
  ##-----------------------------------------------------------
  session.file <-
    base::file(base::paste(output.dirname, 'sessionInfo.txt', sep = '/'),
               open = "wt")
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
