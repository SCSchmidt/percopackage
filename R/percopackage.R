#' percopackage: A package for percolation analysis as a 2D spatial clustering algorithm
#'
#'Three functions are implemented in this package:
#'
#'@section percolate():
#'This function processes this data set to generate the 
#'clusters for a range of percolation radius values. 
#'It also generates and stores a set of data tables for 
#'use by the other functions, as well as optional external 
#'processing. For details see \code{\link{percolate}}
#'
#'@section mapClusters(): 
#'which allows for the mapping of the 
#'clusters together with a shape file as a background. For details see 
#' \code{\link{mapClusters}}
#'
#'@section plotClustFreq(): 
#'which plots three different analyses an png files:
#'
#'a) radius to maximum cluster size, 
#'
#'b) radius to mean cluster size and 
#'
#'c) radius to normalized max. cluster size. 
#'
#'For details see \code{\link{plotClustFreq}}
#'    
#'
#' @docType package
#' @name percopackage
NULL