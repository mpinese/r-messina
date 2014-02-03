#' A class to store the results of a Messina fit.
#'
#' TODO
#' 
#' @name messinaResult-class
#' @seealso \code{\link{messina}}
#' @references TODO
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
NULL


messinaPlot = function(result, i = NULL)
{
	Sample = Value = Class = NULL		# To shut up an R CMD check note for the later use of these in ggplot
	
	if (is.null(i))
	{
		if (!any(result$passed))
		{
			stop("Error: Index of feature to plot was not specified (i == NULL), so best passing feature would be selected.  However, no features passed performance criteria.  Please explicitly specify feature to plot by setting 'i' to a feature index in the plot call.")
		}
		temp.margin = result$margin
		temp.margin[!result$passed] = NA
		i = which.max(temp.margin)
	}
	feature = ifelse(is.null(result$parameters$features), sprintf("index %d", i), result$parameters$features[i])
	
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
		ggtitle(sprintf("Messina Fit: Feature %s", feature)) + coord_flip()
}


calcKaplanMeierEstimates = function(y, x)
{
	xvals = sort(as.character(unique(x)))
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


calcKaplanMeierEstimatesOnBootstrapSample = function(y, x)
{
	samp = sample.int(length(x), replace = TRUE)
	return(calcKaplanMeierEstimates(y[samp,], x[samp]))
}


calcBootstrapKaplanMeierEstimates = function(y, x, nboot)
{
	xvals = sort(as.character(unique(x)))
	results1 = rlply(nboot, calcKaplanMeierEstimatesOnBootstrapSample(y, x))
	results2 = llply(xvals, function(xi) llply(results1, function(ri) ri[[xi]]))
	names(results2) = xvals
	return(results2)
}


calcBootstrapKaplanMeierEstimatesAtTimes = function(y, x, nboot, times = sort(unique(y[,1])))
{
	ests = calcBootstrapKaplanMeierEstimates(y, x, nboot)
	ests_at_times = llply(ests, 
		function(ests_for_x) 
		{
			ests_at_time = laply(ests_for_x, function(est) approx(est$time, est$surv, times, method = "constant", f = 0, rule = 1, ties = "ordered")$y)
			if (nboot == 1)
			{
				ests_at_time = matrix(ests_at_time, nrow = 1)
			}
			result = cbind(times, t(ests_at_time))
			colnames(result) = c("time", paste("surv", 1:nboot, sep = ""))
			result
		}
	)
	names(ests_at_times) = names(ests)
	ests_at_times
}


messinaSurvKMplot = function(y, group, bootstrap_type, bootstrap_ci, nboot)
{
	stopifnot(all(group %in% c(TRUE, FALSE, 0, 1)))
	group = group*1
	
	full_km = calcKaplanMeierEstimates(y, group)

	if (bootstrap_type != "none")
	{
		boot_times = sort(unique(y[,1]))
		boot_km = calcBootstrapKaplanMeierEstimatesAtTimes(y, group, nboot, boot_times)
		names(boot_km) = c("0", "1")

		if (bootstrap_type == "ci")
		{
			boot_km_c = llply(boot_km, function(boots) aaply(boots[,-1,drop=FALSE], 1, median, na.rm = TRUE))
			boot_km_l = llply(boot_km, function(boots) aaply(boots[,-1,drop=FALSE], 1, quantile, probs = 1-bootstrap_ci, na.rm = TRUE))
			boot_km_u = llply(boot_km, function(boots) aaply(boots[,-1,drop=FALSE], 1, quantile, probs = bootstrap_ci, na.rm = TRUE))
		}
		else
		{
			boot_km_sd = llply(boot_km, function(boots) aaply(boots[,-1,drop=FALSE], 1, sd, na.rm = TRUE))
			boot_km_c = llply(c("0", "1"), function(group) rowMeans(boot_km[[group]][,-1]))
			names(boot_km_c) = c("0", "1")
			boot_km_l = llply(c("0", "1"), function(group) boot_km_c[[group]] - boot_km_sd[[group]])
			boot_km_u = llply(c("0", "1"), function(group) boot_km_c[[group]] + boot_km_sd[[group]])
			names(boot_km_l) = c("0", "1")
			names(boot_km_u) = c("0", "1")
		}
		
		boot_data = data.frame(	Time = rep(boot_times, 2), 
								Survival = c(boot_km_c[["0"]], boot_km_c[["1"]]),
								SurvMin = c(boot_km_l[["0"]], boot_km_l[["1"]]),
								SurvMax = c(boot_km_u[["0"]], boot_km_u[["1"]]),
								Group = rep(c("<= Threshold", "> Threshold"), each = length(boot_times)))

		boot_poly_data = data.frame(	Time = rep(rep(c(boot_times[-1], rev(boot_times)[-1]), each = 2), 2),
										Survival = c(c(rep(boot_km_u[["0"]][-1], each = 2)[-1], rep(rev(boot_km_l[["0"]]), each = 2)[-(length(boot_times)*2)]),
													 c(rep(boot_km_u[["1"]][-1], each = 2)[-1], rep(rev(boot_km_l[["1"]]), each = 2)[-(length(boot_times)*2)])),
										Group = rep(c("<= Threshold", "> Threshold"), each = (length(boot_times)-1)*4))
		boot_poly_data = boot_poly_data[!is.na(boot_poly_data$Survival),]
	}
	
	full_data = data.frame(	Time = c(full_km[["0"]]$time, full_km[["1"]]$time), 
							Survival = c(full_km[["0"]]$surv, full_km[["1"]]$surv), 
							Group = factor(c(rep("<= Threshold", nrow(full_km[["0"]])), rep("> Threshold", nrow(full_km[["1"]])))))
		
	if (bootstrap_type != "none")
	{
		theplot = ggplot(data = full_data, aes(x = Time, y = Survival, group = Group, col = Group)) +
			geom_polygon(data = boot_poly_data, mapping = aes(x = Time, y = Survival, group = Group, fill = Group), alpha = 0.2, linetype = 0) + 
			geom_step(direction = "hv") + 
			geom_step(data = boot_data, direction = "hv", alpha = 0.2)
	}
	else
	{
		theplot = ggplot(data = full_data, aes(x = Time, y = Survival, group = Group, col = Group)) +
			geom_step(direction = "hv")
	}
	theplot	
}


messinaSurvObjPlot = function(fit)
{
	plot_data = data.frame(Objective = fit$obj, Cutoff = fit$cutoffs)
	theplot = ggplot(data = plot_data, mapping = aes(x = Cutoff, y = Objective)) + 
		geom_line(alpha = 0.5) + 
		geom_point() + 
		geom_hline(yintercept = fit$obj_min, lty = "dotted") + 
		geom_vline(xintercept = c(fit$threshold, fit$threshold - fit$margin/2, fit$threshold + fit$margin/2), lty = c("solid", "dashed", "dashed"), alpha = c(1, 0.5, 0.5))
	theplot
}


messinaSurvPlot = function(result, i = NULL, bootstrap_type = "none", bootstrap_ci = 0.90, nboot = ifelse(bootstrap_type == "ci", 50/(1-bootstrap_ci), 50))
{
	stopifnot(bootstrap_type %in% c("none", "ci", "stdev"))
	bootstrap_ci = max(bootstrap_ci, 1 - bootstrap_ci)
	stopifnot(bootstrap_ci > 0.5 && bootstrap_ci < 1)
	
	if (is.null(i))
	{
		if (!any(result$passed))
		{
			stop("Error: Index of feature to plot was not specified (i == NULL), so best passing feature would be selected.  However, no features passed performance criteria.  Please explicitly specify feature to plot by setting 'i' to a feature index in the plot call.")
		}
		temp.margin = result$margin
		temp.margin[!result$passed] = NA
		i = which.max(temp.margin)
	}

	x = result$parameters$x[i,]
	y = result$parameters$y
	fit = result$fits[[i]]
	feature = ifelse(is.null(result$parameters$features), sprintf("index %d", i), result$parameters$features[i])
	
	threshold = result$classifier$threshold[i]
	margin = result$margin[i]

	km_plot_threshold = messinaSurvKMplot(y, (x > threshold)*1, bootstrap_type, bootstrap_ci, nboot) + ggtitle("Separation at Threshold")
	km_plot_lower_margin = messinaSurvKMplot(y, (x > threshold - margin/2)*1, bootstrap_type, bootstrap_ci, nboot) + ggtitle("Separation at Lower Margin")
	km_plot_upper_margin = messinaSurvKMplot(y, (x > threshold + margin/2)*1, bootstrap_type, bootstrap_ci, nboot) + ggtitle("Separation at Upper Margin")
	
	obj_plot = messinaSurvObjPlot(fit) + ggtitle("Objective Function")

	grid.arrange(obj_plot, km_plot_threshold, km_plot_lower_margin, km_plot_upper_margin, main = sprintf("MessinaSurv Fit: Feature %s", feature))
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
