#' A class to store the results of a Messina fit.
#'
#' TODO
#' 
#' @name messinaResult-class
#' @seealso \code{\link{messina}}
#' @references TODO
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
NULL


#' Plot the results of a Messina fit.
#'
#' TODO
#' 
#' @param result the result of a Messina fit, either returned by messina or
#'   messinasurv.
#' @param i an optional index of which feature (as ordered by the rows of the
#'   x matrix supplied to messina) to plot.  If missing, the feature with the
#'   largest fit margin is used.
#' @param type the plotting engine to use, either "base" for base graphics, or
#'   "ggplot2" for ggplot2 (the ggplot2 package must be available to use this 
#'   engine).
#' @export
#' @seealso \code{\link{MessinaResult-class}}
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{messinaSurv}}
#' @seealso \code{\link{ggplot2-package}}
#' @references TODO
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
plot.MessinaResult = function(result, ...)
{
	if (missing(i))	i <- NULL
	if (missing(type)) type <- "ggplot2"
	
	Sample = Value = Class = NULL		# To shut up an R CMD check note for the later use of these in ggplot
	
	if (!(type %in% c("base", "ggplot2")))
	{
		stop(sprintf("Error: supplied plot type \"%s\" not recognised.  type must be either \"base\" or \"ggplot2\".", type))
	}
	
	if (type == "ggplot2")
	{
		if (require(ggplot2) == FALSE)
		{
			warning("Warning: ggplot2 library must be available if type = \"ggplot2\".  Falling back to type = \"base\".")
			type = "base"
		}
	}
	
	if (is.null(i))
	{
		temp.margin = result$margin
		temp.margin[!result$passed] = NA
		i = which.max(temp.margin)
	}
	this_x = result$parameters$x[i,]
	this_order = order(this_x)
	this_x = this_x[this_order]
	this_y = result$parameters$y[this_order]
	this_threshold = result$classifier$threshold[i]
	this_margin = result$margin[i]
	ymax = max(c(this_x, this_threshold + this_margin/2))

	this_data = data.frame(Sample = names(this_x), Value = this_x, Class = ordered(this_y*1))

	if (type == "ggplot2")
	{
		#~ ggplot(this_data, aes(x = Value, ymin = 0, fill = Class, colour = Class)) +
			#~ geom_dotplot() +
			#~ geom_vline(xintercept = c(this_threshold - this_margin/2, this_threshold, this_threshold + this_margin/2), linetype = c("dashed", "solid", "dashed"), lwd = 0.7) + 
			#~ ggtitle(result$parameters$features[i])

		ggplot(this_data, aes(x = reorder(Sample, Value), y = Value, fill = Class, colour = Class)) +
			geom_point(stat = "identity") +
			geom_hline(yintercept = c(this_threshold - this_margin/2, this_threshold, this_threshold + this_margin/2), linetype = c("dashed", "solid", "dashed"), lwd = 0.7) + 
			xlab("Sample") +
			ggtitle(result$parameters$features[i]) + coord_flip()
	}
	else
	{		
		barplot(this_x, col = c("green", "red")[this_y + 1], ylim = c(0, ymax))
		abline(h = c(this_threshold - this_margin/2, this_threshold, this_threshold + this_margin/2), lwd = 2, lty = c("dotted", "solid", "dotted"))
	}
}


#' Print a summary of the results of a Messina fit.
#'
#' TODO
#' 
#' @param result an object of class MessinaResult containing the result of a 
#'   messina fit, as returned by the function messina.
#' @export
#' @seealso \code{\link{MessinaResult-class}}
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{plot.MessinaResult}}
#' @references TODO
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
summary.MessinaResult = function(result, ...)
{
	# TODO
}


#' Print the results of a Messina fit, ordered by decreasing margin size.
#'
#' TODO
#' 
#' @param result an object of class MessinaResult containing the result of a 
#'   messina fit, as returned by the function messina.
#' @param maxn the maximum number of top features to display.  If missing, all
#'   features are printed.
#' @param minmar only features with a fit margin of at least minmar will be printed.
#' @export
#' @seealso \code{\link{MessinaResult-class}}
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{plot.MessinaResult}}
#' @references TODO
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
toptable.MessinaResult = function(result, ...)
{
	if (!is.null(maxn))
	{
		ranks = rank(-result$margin, ties.method = "min")
		sel.n = ranks <= maxn
	}
	else	{ sel.n = TRUE }
	
	if (!is.null(minmar))
	{
		sel.mar = result$margin >= minmar
		sel.mar[is.na(sel.mar)] = FALSE
	}
	else	{ sel.mar = TRUE }
	
	sel = sel.n & sel.mar & !is.na(result$margin)
	
	summary = data.frame(	Feature = result$parameters$features[sel],
							Cutoff = result$classifier$threshold[sel],
							HighGroup = result$classifier$posk[sel],
							Margin = result$margin[sel],
							Sensitivity = result$perf$mean[sel,"Sensitivity"],
							Specificity = result$perf$mean[sel,"Specificity"],
							Passed = result$passed[sel])
	summary = summary[order(-summary$Passed, -summary$Margin),]
	rownames(summary) = NULL
	print(summary)
	
	# TODO: add a list of dropped features?  (particularly for margin == NA)
	
	return(invisible(summary))
}
