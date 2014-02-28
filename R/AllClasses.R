setOldClass("Surv")

setClassUnion("MessinaY", c("Surv", "vector"))

#' A class to store the parameters supplied to a messina or messinaSurv analysis.
#'
#' @slot x a matrix of expression values supplied to the messina or messinaSurv
#'   functions.  Features are in rows, samples in columns.
#' @slot y either a vector of class membership indicators (0/1 or TRUE/FALSE), for
#'   the messina case, or a Surv object for the messinaSurv case.  In either case,
#'   the length of y is equal to the number of columns of x, and each entry of y
#'   should match the corresponding sample column of x.
#' @slot features a character vector of feature ids, matching the rows of x.
#' @slot samples a character vector of sample ids, matching the columns of x and 
#'   entries of y.
#' @slot perf_requirement a list of performance requirements.  For messina results,
#'   contains named entries "min_sensitivity" and "min_specificity".  For messinaSurv
#'   results, contains named entries "objective_type" and "min_objective".
#' @slot minimum_group_fraction the size, relative to the full sample size, of the
#'   smallest subgroup that may be defined by a threshold.
#' @slot training_fraction the fraction of samples used for training in each bootstrap
#'   round.
#' @slot num_bootstraps the number of bootstrap iterations to perform.
#' @slot prng_seed the PRNG seed used to initialize the PRNG before analysis.
#' 
#' @aliases MessinaParameters-class
#'
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{messinaSurv}}
#' @seealso \code{\link{MessinaResult-class}}
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
.MessinaParameters <- setClass(	"MessinaParameters",
								slots = c(	x = "matrix",
											y = "MessinaY",
											features = "vector",
											samples = "vector",
											perf_requirement = "list",		# For MessinaSurvResult, list with names objective_type, min_objective.  For MessinaClassResult, list with names min_sensitivity, min_specificity
											minimum_group_fraction = "numeric",
											training_fraction = "numeric",
											num_bootstraps = "integer",
											prng_seed = "integer"))
											
#' A class to store the individual messina or messinaSurv fits to a
#' dataset.
#'
#' @slot summary a data frame containing summary performance measures
#'   for each feature, features in rows.  The user-relevant columns are:
#'     \describe{
#'       \item{passed}{did this feature pass the user requirements?  A boolean.}
#'       \item{type}{the type of classifier that was fit}
#'       \item{threshold}{the threshold expression value of the classifier}
#'       \item{posk}{the direction of the classifier}
#'       \item{margin}{the expression margin of the classifier}
#'     }
#' @slot objective_surfaces a list of length equal to the number of features.
#'   each list entry contains a data frame of the objective function values
#'   at each threshold (cutoff) tested.  Currently only populated for 
#'   messinaSurv fits.
#' 
#' @aliases MessinaFits-class
#'
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{messinaSurv}}
#' @seealso \code{\link{MessinaResult-class}}
#' @seealso \code{\link{MessinaParameters-class}}
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
.MessinaFits <- setClass("MessinaFits",
						slots = c(	summary = "data.frame",					# columns: passed, type, threshold, posk, margin, ptrue, psuccessful
									objective_surfaces = "list"))			# list of data.frames, each with columns: cutoff, objective (MessinaSurvResult), or cutoff, sensitivity, specificity (MessinaClassResult)

#' A class to store the results of a messina or messinaSurv analysis.
#'
#' @slot problem_type The variant of the messina algorithm used, either
#'   "classification" for the classification case (fit using the function
#'   messina, or "survival" for the outcome case (fit using the function
#'   messinaSurv).
#' @slot parameters An object of class MessinaParameters, containing
#'   input data and parameters for the algorithm.
#' @slot perf_estimates A data frame of summary performance estimates
#'   (evaluated on many out-of-bag sample draws), with one row per feature in the
#'   data matrix supplied to the fit functions (either messina or 
#'   messinaSurv).  For a messina fit, this contains 10 columns:
#'   Mean TPR, Mean FPR, Mean TNR, Mean FNR, Variance of TPR, Variance
#'   of FPR, Variance of TNR, Variance of FNR, Mean sensitivity, Mean
#'   specificity.  For a messinaSurv fit, this contains a single column,
#'   of the mean objective value for that row's feature.
#' @slot fits An object of class MessinaFits, containing details of the
#'   fits for each feature.
#'
#' @aliases MessinaResult-class MessinaClassResult-class MessinaSurvResult-class
#'
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{messinaSurv}}
#' @seealso \code{\link{MessinaParameters-class}}
#' @seealso \code{\link{MessinaFits-class}}
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
.MessinaResult <- setClass(	"MessinaResult", 
							slots = c(	problem_type = "character",
										parameters = "MessinaParameters",
										perf_estimates = "data.frame",		# Columns for MessinaClassResult: mean_tpr, mean_fpr, mean_tnr, mean_fnr, var_tpr, var_fpr, var_tnr, var_fnr, mean_sens, mean_spec
																			# Columns for MessinaSurvResult: mean_obj
										fits = "MessinaFits"))

#' @export
.MessinaClassResult <- setClass("MessinaClassResult", contains = "MessinaResult")

#' @export
.MessinaSurvResult <- setClass(	"MessinaSurvResult", contains = "MessinaResult")
