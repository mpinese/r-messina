setGeneric("topTable", function(object, ...) standardGeneric("topTable"))
#if (!isGeneric("topTable"))		{ setGeneric("topTable", function(object, ...) standardGeneric("topTable")) }
# TODO: The above is safer, but currently doesn't work with roxygen2 (I assume it skips this code).
# Figure out what to do.

#' @export
setMethod("topTable", signature = signature(object = "MessinaResult"), definition = function(object, ...) messinaTopTable(object@fits, ...))


messinaTopTable = function(fits, n = 10)
{
	summary = fits@summary
	n = min(n, nrow(summary))
	summary = summary[order(-summary$margin),]
	summary = summary[,c("passed", "type", "threshold", "posk", "margin")]
	summary$posk = c(-1, 1)[summary$posk + 1]
	colnames(summary) = c("Passed Requirements", "Classifier Type", "Threshold Value", "Direction", "Margin")
	print(summary[1:n,])
	invisible(summary)
}
