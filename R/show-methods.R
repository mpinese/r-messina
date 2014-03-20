# show-methods: Methods for the show generic on Messina objects.
# 
# Copyright 2014 Mark Pinese
#
# This file is distributed under the terms of the Eclipse Public 
# License v1.0, available at:
# https://www.eclipse.org/org/documents/epl-v10.html


#' Generic show methods for Messina objects.
#'
#' Generic show methods for Messina objects.
#' 
#' For details of the objects and their generation, see the relevant class documentation,
#' and entries for the main functions \code{\link{messina}}, \code{\link{messinaDE}}, and 
#' \code{\link{messinaSurv}}, 
#' 
#' @importFrom methods show
#' 
#' @export
#' @seealso \code{\link{MessinaResult-class}}
#' @seealso \code{\link{MessinaParameters-class}}
#' @seealso \code{\link{MessinaFits-class}}
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{messinaDE}}
#' @seealso \code{\link{messinaSurv}}
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
#' 
#' @importFrom methods show
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
#' 
#' @importFrom methods show
setMethod("show", signature = "MessinaFits", definition = function(object) 
{
	cat("  An object of class MessinaFits\n", sep = "")	
	cat("  ", sum(object@summary$passed), " / ", length(object@summary$passed), sprintf(" features passed performance requirements (%.2f%%)\n", mean(object@summary$passed)*100), sep = "")
	cat("  Top features:\n")
	
	messinaTopResults(object, 10)
	
	invisible(NULL)
})
