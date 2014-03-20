# messinTopResults.R: Functions to display a summary table of Messina results
# 
# Copyright 2014 Mark Pinese
#
# This file is distributed under the terms of the Eclipse Public 
# License v1.0, available at:
# https://www.eclipse.org/org/documents/epl-v10.html


#' Display a summary of the top results from a Messina analysis
#' 
#' Sorts the summary results of a Messina analysis in decreasing order of classifier margin,
#' and displays the top n.  The full sorted data.frame is invisibly returned.
#'
#' The displayed data.frame has the following columns.  Users are encouraged to consult the vignette
#' for a tutorial on how to interpret these results for classification and gene expression tasks.
#' \describe{
#'   \item{"Rowname"}{The feature ID.  If the data x supplied to the messina function was an ExpressionSet
#'     object, the featureName of the relevant feature of x.  Otherwise, if x was a matrix with row names, 
#'     the row name of the corresponding entry in the matrix, or if x was a matrix without row names, F<n>, 
#'     where <n> is the row number of the corresponding row of x.}
#'   \item{"Passed Requirements"}{Logical: when its performance was assessed by bootstrapping, did this feature 
#'     pass the user-supplied performance requirements on out-of-bag data?}
#'   \item{"Classifier Type"}{A string indicating the type of single-gene classifier that Messina fit to 
#'     this feature.  Valid values are given below, but for most users only the Threshold type is relevant,
#'     with the others being only of diagnostic relevance.
#'     \describe{
#'       \item{"Threshold"}{A threshold classifer: samples with feature signal at or below the threshold are in one
#'         group; samples with feature signal above the threshold are in the other.  This is the main result
#'         of interest in a Messina analysis, and other classifier types are more of diagnostic interest.}
#'       \item{"Random"}{A random (also known as Zero-Rule) classifier.  In this case, the feature did not contain
#'         sufficient information to construct a good classifier, but the performance requirements were so
#'         lenient that simple guessing of an unknown sample's class based on marginal probabilities was enough
#'         to satisfy them.  The presence of these 'fits' in the top results is indicative of too lenient 
#'         performance requirements, or a dataset with no predictive value for the classes of interest (at least
#'         for single-feature threshold classifiers).}
#'       \item{"OneClass"}{All samples are always called as a single class, and this strategy is sufficient to
#'         satisfy the supplied performance requirements.  Similar to the "Random" type, the presence of these
#'         results are an indicator of too lenient performance requirements.}
#'       \item{"NA"}{The feature was not successfully fit.  Seen as an indicator of failed fitting in MessinaSurv
#'         analyses only, where the Random and OneClass defaults are not applicable.}
#'     }
#'   }
#'   \item{"Threshold Value"}{For a Threshold classifier, the value of the optimal threshold selected by the
#'     algorithm.  This is the value to use as a cutoff in separating the samples into two classes, either
#'     "Group 0" and "Group 1", or "Long surviors" and "Short survivors".}
#'   \item{"Direction"}{The direction of the threshold classifier.  Can take values of either -1 or 1.  If -1,
#'     samples with expression above the threshold are in group 1 (/TRUE), or have shorter survival times.
#'     If 1, samples with expression value above the threshold are in group 0 (/FALSE), and have longer survival
#'     times.}
#'   \item{"Margin"}{The value of the threshold classifier's margin.  This is the primary measure of fit strength
#'     in a Messina analysis: a higher margin indicates stronger robustness to noise and experimental variations
#'     in a classification context, and a higher likelihood of differential expression in a gene expression
#'     context.}
#' }
#'
#' @param result the result returned by a call to messina, messinaDE, or messinaSurv.
#' @param n the maximum number of top hits to display (default 10).
#'
#' @return (invisible) the full table of hits, as a data.frame sorted in order 
#'   of decreasing margin.
#'
#' @export
#'
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{messinaDE}}
#' @seealso \code{\link{messinaSurv}}
#' @seealso \code{\link{MessinaClassResult-class}}
#'
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
#'
#' @examples
#' ## Load some example data
#' library(antiProfilesData)
#' data(apColonData)
#' 
#' x = exprs(apColonData)
#' y = pData(apColonData)$SubType
#' 
#' ## Subset the data to only tumour and normal samples
#' sel = y %in% c("normal", "tumor")
#' x = x[,sel]
#' y = y[sel]
#' 
#' ## Find differentially-expressed probesets.  Allow a sample misattribution rate of
#' ## at most 20%.
#' fit = messina(x, y == "tumor", min_sens = 0.95, min_spec = 0.85)
#'
#' ## Print the 20 probesets with the strongest evidence for differential expression
#' ## between tumour and normal.  Save the full table of summary results for later use.
#' summary_table = messinaTopResults(fit, 20)
#'
#' ## Access the top five probesets in the table
#' summary_table[1:5,]
#'
#' ## Examine the summary results for particular probes
#' summary_table[c("204719_at", "207502_at"),]
messinaTopResults = function(result, n = 10)
{
	if (is(result, "MessinaResult"))
	{
		fits = result@fits
	}
	else if (is(result, "MessinaFits"))
	{
		fits = result
	}
	else
	{
		stop("messinaTopResults: Result object must be of class MessinaClassResult, MessinaSurvResult, or MessinaFits")
	}

	summary = fits@summary
	n = min(n, nrow(summary))
	summary = summary[order(-summary$passed, -summary$margin),]
	summary = summary[,c("passed", "type", "threshold", "posk", "margin")]
	summary$posk = c(-1, 1)[summary$posk + 1]
	colnames(summary) = c("Passed Requirements", "Classifier Type", "Threshold Value", "Direction", "Margin")
	print(summary[1:n,])
	invisible(summary)
}
