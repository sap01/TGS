#' Create a .sif file from given adjacency matrix
#'
#' Given a network adjacency matrix, creates an equivalent '.sif' file that
#' is readable in Cytoscape. SIF file extension stands for Simple Interaction
#' Format.
#' Ref: http://manual.cytoscape.org/en/3.4.0/Supported_Network_File_Formats.html#sif-format
#' E.g., if the input network adjacency matrix is as follows: (rows = src nodes, cols = tgt nodes,
#' (A, B) = 0 and 1 implies that the edge 'A->B' does not exist and does exist, resp.)
#' src\ tgt  A B C D
#' A         0 1 1 0
#' B         0 0 0 0
#' C         1 0 0 0
#' C         0 0 0 0
#' then the output SIF file will contain the following four lines:
#' A  B C
#' B
#' C  A
#' D
#' The way Cytoscape interprets each line in this SIF file is as follows:
#' <Source node name> <Target node1 name (if any)>  <Target node2 name (if any)> ...
#' Note that the values are delimited by tab.
#'
#' @param adj.mx adjacency matrix that needs to be converted
#' @param output.dirname name of the output directory where the file will be stored
#'
#' @export
adjmxToSif <- function(adj.mx, output.dirname = "./OUTPUT")
{
  if(!base::is.matrix(adj.mx))
  {
    base::stop("Error in adjmxToSif adj.mx is not a matrix")
  }
  ## Open an output file connection in write mode
  output.sif <- base::file(base::paste(output.dirname, 'net.sif', sep = '/'), 'w')

  for (src.node.idx in 1:nrow(adj.mx))
  {
    ## 'pd' stands for Protein-DNA interaction type in Cytoscape.
    line.to.write <- base::paste(base::rownames(adj.mx)[src.node.idx], 'pd', sep = '\t')

    for (tgt.node.idx in 1:ncol(adj.mx))
    {
      if (adj.mx[src.node.idx, tgt.node.idx] == 1)
      {
        line.to.write <- base::paste(line.to.write, base::colnames(adj.mx)[tgt.node.idx], sep = '\t')
      }
    }
    base::rm(tgt.node.idx)

    base::cat(line.to.write, file = output.sif, '\n')
  }
  base::rm(src.node.idx)

  ## Close the output file connection
  base::close(output.sif)
}
