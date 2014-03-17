#' Single-gene classifiers and outlier-resistant detection of differential expression 
#' for two-group and survival problems.
#'
#' Messina is a collection of algorithms for constructing optimally robust
#' single-gene classifiers, and for identifying differential expression in
#' the presence of outliers or unknown sample subgroups.  The methods have
#' application in identifying lead features to develop into clinical 
#' tests (both diagnostic and prognostic), and in identifying differential
#' expression when a fraction of samples show unusual patterns of expression.
#'
#' @useDynLib messina
#' @name messina-package
#' @docType package
#' @title The Messina package for classification and outlier differential expression.
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
#' @references Pinese:2009 Pinese M, Scarlett CJ, Kench JG, et al. (2009)
#'   Messina: A Novel Analysis Tool to Identify Biologically Relevant 
#'   Molecules in Disease.  PLoS ONE 4(4): e5337.  \url{doi:10.1371/journal.pone.0005337}
#' @keywords package
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{messinaDE}}
#' @seealso \code{\link{messinaSurv}}
NULL
