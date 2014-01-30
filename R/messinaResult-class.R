#' A class to store the results of a Messina fit.
#'
#' TODO
#' 
#' @name messinaResult-class
#' @seealso \code{\link{messina}}
#' @references TODO
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
NULL


messinaPlot = function(result, i)
{
	Sample = Value = Class = NULL		# To shut up an R CMD check note for the later use of these in ggplot
	
	if (is.null(i))
	{
		temp.margin = result$margin
		temp.margin[!result$passed] = NA
		i = which.max(temp.margin)
	}
	x = result$parameters$x[i,]
	samples = names(x)
	if (is.null(samples))	samples = seq_along(x)
	order = order(x)
	x = x[order]
	samples = samples[order]
	y = result$parameters$y[order]
	threshold = result$classifier$threshold[i]
	margin = result$margin[i]

	data = data.frame(Sample = samples, Value = x, Class = ordered(y*1))

	ggplot(data, aes(x = reorder(Sample, Value), y = Value, fill = Class, colour = Class)) +
		geom_point(stat = "identity") +
		geom_hline(yintercept = c(threshold - margin/2, threshold, threshold + margin/2), linetype = c("dashed", "solid", "dashed"), lwd = 0.7) + 
		xlab("Sample") +
		ggtitle(result$parameters$features[i]) + coord_flip()
}


getKaplanMeierEstimates = function(y, x)
{
	xvals = as.character(unique(x))
	result = llply(xvals, 
		function(xi)
		{
			fit = survfit(y[x == xi,] ~ 1)
			return(data.frame(time = fit$time, surv = fit$surv))
		}
	)
	names(result) = xvals
	return(result)
}


getKaplanMeierEstimatesOnBootstrapSample = function(y, x)
{
	samp = sample.int(length(x), replace = TRUE)
	return(getKaplanMeierEstimates(y[samp,], x[samp]))
}


getBootstrapKaplanMeierEstimates = function(y, x, nboot)
{
	xvals = as.character(unique(x))
	results1 = rlply(nboot, getKaplanMeierEstimatesOnBootstrapSample(y, x))
	results2 = llply(xvals, function(xi) ldply(results1, function(ri) ri[[xi]]))
	return(results2)
}

getKaplanMeierEstimates(data.ys, data.x[1,] > 15)
getKaplanMeierEstimatesOnBootstrapSample(data.ys, data.x[1,] > 15)
debug(getBootstrapKaplanMeierEstimates)
getBootstrapKaplanMeierEstimates(data.ys, data.x[1,] > 15, nboot = 2)
plot(getBootstrapKaplanMeierEstimates(data.ys, data.x[1,] > 15, nboot = 100)[[1]], type = "s")

messinaSurvPlot = function(result, i, nboot = 0)
{
	if (is.null(i))
	{
		temp.margin = result$margin
		temp.margin[!result$passed] = NA
		i = which.max(temp.margin)
	}

	x = result$parameters$x[i,]
	y = result$parameters$y
	
	threshold = result$classifier$threshold[i]
	margin = result$margin[i]
	posk = result$classifier$posk[i]

	full_km_lower = getKaplanMeierEstimates(y, xor(x > (threshold - margin/2), !posk)*1)
	full_km_upper = getKaplanMeierEstimates(y, xor(x > (threshold + margin/2), !posk)*1)
	full_km_threshold = getKaplanMeierEstimates(y, xor(x > threshold, !posk)*1)

	boot_km_lower = getBootstrapKaplanMeierEstimates(y, xor(x > (threshold - margin/2), !posk)*1, nboot)
	full_km_upper = getBootstrapKaplanMeierEstimates(y, xor(x > (threshold + margin/2), !posk)*1, nboot)
	full_km_threshold = getBootstrapKaplanMeierEstimates(y, xor(x > threshold, !posk)*1, nboot)

	# "s"
	
	# Plots: Objective, KM at optimal threshold,
	# KM at lower limit, KM at upper limit.
}


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
	if (result$problem == "classification")
	{
		messinaPlot(result, ...)
	}
	else if (result$problem == "survival")
	{
		messinaSurvPlot(result, ...)
	}
	else
	{
		stop(sprintf("Unknown Messina problem type: \"%s\"", result$problem))
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
toptable.MessinaResult = function(result, maxn, minmar)
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












	
#~ plotObj = function(x, y, func, inv = FALSE, ...)
#~ {
	#~ plotx = seq(min(x), max(x), length.out = 500)
	#~ if (inv == TRUE)
	#~ {
		#~ ploty = sapply(plotx, function(xi) messinaSurvObjectiveFunc(x < xi, y, func))
	#~ }
	#~ else
	#~ {
		#~ ploty = sapply(plotx, function(xi) messinaSurvObjectiveFunc(x > xi, y, func))
	#~ }
	
	#~ plot(ploty ~ plotx, ...)
#~ }


#~ plotTrain = function(x, y, fit, ...)
#~ {
	#~ if (is.na(fit$threshold))
	#~ {
		#~ plot.new()
		#~ plot.new()
		#~ plot.new()
	#~ }
	#~ else
	#~ {
		#~ xt = x * c(-1, 1)[1*(fit$direction == 1) + 1]
		#~ tt = fit$threshold * c(-1, 1)[1*(fit$direction == 1) + 1]
		#~ plot(survfit(y ~ I(xt > tt)), main = "Separation at Threshold", xlab = "Time", ylab = "Surviving fraction", col = c("red", "blue"))
		#~ tt = (fit$threshold - fit$margin/2) * c(-1, 1)[1*(fit$direction == 1) + 1]
		#~ plot(survfit(y ~ I(xt > tt)), main = "Separation at Lower Bound", xlab = "Time", ylab = "Surviving fraction", col = c("red", "blue"))
		#~ tt = (fit$threshold + fit$margin/2) * c(-1, 1)[1*(fit$direction == 1) + 1]
		#~ plot(survfit(y ~ I(xt > tt)), main = "Separation at Upper Bound", xlab = "Time", ylab = "Surviving fraction", col = c("red", "blue"))
	#~ }
	#~ plot(fit$obj ~ fit$cutoffs, main = "Objective Surface", xlab = "Threshold", ylab = "Objective", type = "o", ...)
	#~ abline(v = c(fit$threshold, fit$threshold - fit$margin/2, fit$threshold + fit$margin/2), lty = c("solid", "dashed", "dashed"))
	#~ abline(h = fit$obj.min, lty = "dotted", col = "grey")
#~ }


#~ summarizeResult = function(result)
#~ {
	#~ temp = result$summary[result$summary$PassedObj == TRUE,]
	#~ temp = temp[order(-temp$Margin),]
	#~ temp
#~ }

