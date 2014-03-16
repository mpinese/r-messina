#' TODO
#' 
#' TODO
#' 
#' @export
#' @seealso \code{\link{MessinaResult-class}}
#' @seealso \code{\link{messinaSurv}}
#' @seealso \code{\link{messina}}
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
#' @rdname show-methods
setMethod("show", signature = "MessinaResult", definition = function(object) 
{
	cat("An object of class ", class(object), "\n\n", sep = "")
	cat("Problem type:", object@problem_type, "\n", sep = "")
	cat("Parameters:\n")
	show(object@parameters)
	cat("\n")
	cat("Summary of results:\n")
	show(object@fits)
	invisible(NULL)
})


#' @export
#' @rdname show-methods
setMethod("show", signature = "MessinaParameters", definition = function(object) 
{
	cat("  An object of class MessinaParameters\n", sep = "")
	cat("  ", nrow(object@x), " features, ", ncol(object@x), " samples.\n", sep = "")
	if ("min_sensitivity" %in% names(object@perf_requirement))
	{
		cat("  Objective type: sensitivity/specificity.  Minimum sensitivity: ", object@perf_requirement$min_sensitivity, "  Minimum specificity: ", object@perf_requirement$min_specificity, "\n", sep = "")
	}
	else
	{
		cat("  Objective type: survival (", object@perf_requirement$objective_type, ").  Minimum objective value: ", object@perf_requirement$min_objective, "\n", sep = "")
	}
	cat("  Minimum group fraction: ", object@minimum_group_fraction, "\n", sep = "")
	cat("  Training fraction: ", object@training_fraction, "\n", sep = "")
	cat("  Number of bootstraps: ", object@num_bootstraps, "\n", sep = "")
	cat("  Random seed: ", object@prng_seed, "\n", sep = "")
	
	invisible(NULL)
})


#' @export
#' @rdname show-methods
setMethod("show", signature = "MessinaFits", definition = function(object) 
{
	cat("  An object of class MessinaFits\n", sep = "")	
	cat("  ", sum(object@summary$passed), " / ", length(object@summary$passed), sprintf(" features passed performance requirements (%.2f%%)\n", mean(object@summary$passed)*100), sep = "")
	cat("  Top features:\n")
	
	messinaTopResults(object, 10)
	
	invisible(NULL)
})
