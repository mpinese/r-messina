# plot-methods.R: Methods for the plot generic on MessinaResult objects
# 
# Copyright 2014 Mark Pinese
#
# This file is distributed under the terms of the Eclipse Public 
# License v1.0, available at:
# https://www.eclipse.org/org/documents/epl-v10.html


#' Plot the results of a Messina analysis on a classification / differential expression problem.
#' 
#' Produces a separate plot for each supplied feature index (either as an index into the expression
#' data x as-supplied, or as an index into the features sorted by Messina margin, depending on the
#' value of sort_features), showing sample expression levels, group membership, threshold value,
#' and margin locations.  Two different types of plots can be produced.  See the vignette for
#' examples.
#'
# @usage plot(object, indices = c(1), sort_features = TRUE, plot_type = "bar", ...)
#'
#' @param object the result of a Messina analysis, as returned by functions \code{\link{messina}}
#'   or \code{\link{messinaDE}}.
#' @param indices a vector of indices of features to plot.  If sort_features == FALSE, the indices
#'   are into the unsorted features, as originally supplied in x supplied to messina or messinaDE.
#'   If sort_features == TRUE, features are first sorted in order of decreasing margin, and then 
#'   the indices in this parameter are plotted.  For example, if indices == 2 and sort_features == FALSE,
#'   the second feature in x will be plotted.  However, if sort_features == TRUE, the feature with
#'   the second best classifier margin will be plotted.
#' @param sort_features a boolean indicating whether to sort features by decreasing margin size
#'   before selecting from indices.  This affects the interpretation of the parameter 'indices'; for
#'   more details see the description of that parameter.
#' @param plot_type a string giving the type of plot to produce, either "point" or "bar".  "bar"
#'   is the default, and shows expression levels as horizontal bars.  Although this representation
#'   is familiar, it can be misleading in the case of log-transformed data.  In that case, the 
#'   "point" plot type is preferable.
#' 
#' @aliases plot,MessinaClassResult-method
#' @aliases plot,MessinaClassResult,missing-method
#'
#' @export
#' @seealso \code{\link{MessinaClassResult-class}}
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{messinaDE}}
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
#' ## Run Messina to rank probesets on their classification ability, with
#' ## classifiers needing to meet a minimum sensitivity of 0.95, and minimum
#' ## specificity of 0.85.
#' fit = messina(x, y == "tumor", min_sens = 0.95, min_spec = 0.85)
#'
#' ## Make bar plots of the five best fits
#' plot(fit, indices = 1:5, sort_features = TRUE, plot_type = "bar")
#'
#' ## Make a point plot of the fit to the 10th feature
#' plot(fit, indices = 10, sort_features = FALSE, plot_type = "point")
setMethod("plot", signature = signature(x = "MessinaClassResult", y = "missing"), definition = function(x, y, ...) messinaClassPlot(object = x, ...))



#' Plot the results of a Messina analysis on a survival problem.
#'
#' Plots diagnostic and performance information for fits in a MessinaSurvResult object, as returned by
#' \code{\link{messinaSurv}}.
#'
#' For each feature index given by indices, produces four plots:
#' \describe{
#'   \item{"Objective function"}{A plot of the value of the objective function over all possible thresholds.
#'     Each sample is represented by a point on the objective function trace.  The selected threshold, if
#'     any, is shown by a solid vertical line, and the margins by dotted vertical lines on either side of
#'     it.  The minimum values of the objective function specified by the user are shown as horizontal
#'     dotted lines.  This plot is useful for assessing fit stability, particularly for the "coxcoef" and
#'     "reltau" objective functions, which can be unstable at low or high threshold values.  See \code{\link{messinaSurv}}
#'     for details.}
#'   \item{"Separation performance at threshold"}{This Kaplan-Meier plot shows two traces, showing the 
#'     outcomes of the two subgroups in the cohort defined by whether the plotted feature is above or
#'     below the threshold.  Optionally (if bootstrap_type != "none"), the KM traces will be surrounded
#'     by shaded regions that represent either +/- 1 SD (bootstrap_type == "stdev") or a bootstrap_ci
#'     confidence interval (bootstrap_type == "ci").}
#'   \item{"Separation performance at lower margin"}{This plot is identical to the above, except that the
#'     performance when the lower margin is used to separate the sample groups is shown.}
#'   \item{"Separation performance at lower margin"}{This plot is identical to the above, except that the
#'     performance when the upper margin is used to separate the sample groups is shown.  These last 
#'     two plots give an indication of the robustness of the MessinaSurv fit at its extremes.}}
#'
#' The Kaplan-Meier plots may optionally display bootstrap bands, if bootstrap_type != "none".  Note that
#' the calculation of bootstrap bands is computationally-intensive, and this function will by default use
#' multiprocessing to speed calculations if doMC is loaded and more than one core registered for use.
#' For examples of the plots and their interpretation, see the vignette.
#'
# @usage plot(object, indices = c(1), sort_features = TRUE, bootstrap_type = "none", bootstrap_ci = 0.90, nboot = ifelse(bootstrap_type == "ci", 50/(1-bootstrap_ci), 50), parallel = ("doMC" %in% .packages()) && (getDoParWorkers() > 1)), ...)
#'
# @inheritParams plot,MessinaClassResult,missing-method
#' @param object the result of a Messina analysis, as returned by functions \code{\link{messina}}
#'   or \code{\link{messinaDE}}.
#' @param indices a vector of indices of features to plot.  If sort_features == FALSE, the indices
#'   are into the unsorted features, as originally supplied in x supplied to messina or messinaDE.
#'   If sort_features == TRUE, features are first sorted in order of decreasing margin, and then 
#'   the indices in this parameter are plotted.  For example, if indices == 2 and sort_features == FALSE,
#'   the second feature in x will be plotted.  However, if sort_features == TRUE, the feature with
#'   the second best classifier margin will be plotted.
#' @param sort_features a boolean indicating whether to sort features by decreasing margin size
#'   before selecting from indices.  This affects the interpretation of the parameter 'indices'; for
#'   more details see the description of that parameter.
#' @param bootstrap_type a string giving the type of bootstrap error band to produce on the survival prediction
#'   plots.  Can take three values: "none", "stdev", and "ci".  "none", the default, plots no error bands.
#'   "stdev" performs multiple rounds of Kaplan-Meier curve estimation on bootstrap samples,
#'   and plots prediction bands corresponding to +/- 1 bootstrap standard deviation from the mean.  "ci"
#'   performs bootstrapping as per "stdev", and plots prediction bands corresponding to the bootstrap_ci
#'   intervals.
#' @param bootstrap_ci a value in (0.5, 1) giving the confidence interval for bootstrap_type == "ci".  
#'   Ignored otherwise.  Default 0.9 for 90\% confidence intervals.
#' @param nboot the number of bootstrap iterations to perform for calculations.  Set to a reasonable default
#'   taking into account bootstrap_type and bootstrap_ci, so ordinarily does not need to be specified by
#'   the user.
#' @param parallel a logical indicating whether multiprocessing using doMC should be used for the bootstrap
#'   calculations.  If NULL, multiprocessing will be used if doMC is loaded and more than one parallel
#'   worker is registered.
#' 
#' @aliases plot,MessinaSurvResult-method
#' @aliases plot,MessinaSurvResult,missing-method
#'
#' @export
#' @seealso \code{\link{MessinaSurvResult-class}}
#' @seealso \code{\link{messinaSurv}}
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
#'
#' @examples
#' \dontrun{
#' ## Load a subset of the TCGA renal clear cell carcinoma data
#' ## as an example.
#' data(tcga_kirc_example)
#' 
#' ## Run the messinaSurv analysis on these data.  Use a tau
#' ## objective, with a minimum performance of 0.6.  Note that
#' ## messinaSurv analyses are very computationally-intensive,
#' ## so multicore use with doMC loaded and parallel = TRUE is
#' ## strongly recommended.
#' library(doMC)
#' registerDoMC(32)
#' fit = messinaSurv(kirc.exprs, kirc.surv, obj_func = "tau", obj_min = 0.6, parallel = TRUE)
#'
#' ## Plot the three best features found by Messina
#' plot(fit, indices = 1:3)
#'
#' ## Plot the best feature found by Messina, with 90% confidence bands.
#' ## Note that the bootstrap iterations can be slow, so it is 
#' ## recommended that multiple cores are used, with doMC loaded 
#' ## and parallel = TRUE
#' plot(fit, indices = 1, bootstrap_type = "ci", bootstrap_ci = 0.9, parallel = TRUE)
#'
#' ## Plot the Messina fit of the 10th feature in the dataset, with
#' ## +/- 1 standard deviation bands.
#' plot(fit, indices = 10, sort_features = FALSE, bootstrap_type = "stdev")
#' }
setMethod("plot", signature = signature(x = "MessinaSurvResult", y = "missing"), definition = function(x, y, ...) messinaSurvPlot(object = x, ...))


#' @import ggplot2
#' @importFrom grid grid.newpage
messinaClassPlot = function(object, indices = c(1), sort_features = TRUE, plot_type = "bar")
{
	Sample = Value = Class = NULL		# To shut up an R CMD check note for the later use of these in ggplot
	
	n = nrow(object@parameters@x)
	
	if (any(indices > n) || any(indices < 1))
	{
		warning("Warning: Some feature indices were out of range.  Dropping invalid indices.")
		indices = indices[indices >= 1 & indices < n]
	}
	
	if (length(indices) == 0)
	{
		return()
	}
	
	fit_summary = object@fits@summary
	
	feature_perm = 1:n
	if (sort_features)
	{
		feature_perm = c((1:n)[fit_summary$passed][order(-fit_summary$margin[fit_summary$passed])], (1:n)[!fit_summary$passed])
	}
	
	indices = feature_perm[indices]

	for (i in indices)
	{
		grid.newpage()

		feature = object@parameters@features[i]
		x = object@parameters@x[i,]

		order = order(x)
		x = x[order]
		samples = object@parameters@samples[order]
		y = object@parameters@y[order]
		threshold = fit_summary$threshold[i]
		margin = fit_summary$margin[i]

		data = data.frame(Sample = samples, Value = x, Class = ordered(y*1))

		theplot = ggplot(data, aes(x = reorder(Sample, Value), y = Value, fill = Class, colour = Class)) +
			xlab("Sample") +
			ggtitle(sprintf("MessinaClass Fit: Feature %s", feature)) + 
			coord_flip()

		if (plot_type == "point")		{ theplot = theplot + geom_point(stat = "identity") }
		else if (plot_type == "bar")	{ theplot = theplot + geom_bar(stat = "identity", width = 0.5) }
		else 							{ stop(sprintf("Invalid Messina plot type: \"%s\"", plot_type)) }	

		if (!is.na(threshold))
		{
			theplot = theplot + 
				geom_hline(yintercept = threshold, linetype = "solid", lwd = 0.7)
				
			if (!is.na(margin) && margin != 0)
			{
				theplot = theplot + 
					geom_hline(yintercept = c(threshold - margin/2, threshold + margin/2), linetype = "dashed", lwd = 0.7)
			}
		}
		
		print(theplot)
	}
}


#' @import ggplot2
#' @importFrom grid grid.newpage viewport pushViewport popViewport grid.layout
messinaSurvPlot = function(object, indices = c(1), sort_features = TRUE, bootstrap_type = "none", bootstrap_ci = 0.90, nboot = ifelse(bootstrap_type == "ci", 50/(1-bootstrap_ci), 50), parallel = NULL)
{
	if (!(bootstrap_type %in% c("none", "ci", "stdev")))
	{
		stop(sprintf("Error: Invalid bootstrap type, \"%s\".  Valid values for bootstrap_type are \"none\", \"ci\", and \"stdev\".", bootstrap_type))
	}

	if (bootstrap_ci <= 0 || bootstrap_ci >= 1)
	{
		stop(sprintf("Error: Invalid value for bootstrap confidence interval percentile, \"%s\".  Valid values for bootstrap_ci are in (0, 1)  (and usually close to 1!).", str(bootstrap_ci)))
	}
	
	if (bootstrap_ci < 0.7)
	{
		warning(sprintf("Warning: Unusually low value for bootstrap_ci, %f.  Are you sure you want %.0f%% confidence intervals?  Conventionally, bootstrap_ci should be at least 0.9.\nContinuing anyway.", bootstrap_ci, bootstrap_ci*100))
	}

	if (is.null(parallel))
	{
		if ("doMC" %in% .packages())
		{
			parallel = (getDoParWorkers() > 1)
		}
		else
		{
			parallel = FALSE
		}
	}

	nboot = as.integer(round(nboot))

	n = nrow(object@parameters@x)
	
	if (any(indices > n) || any(indices < 1))
	{
		warning("Warning: Some feature indices were out of range.  Dropping invalid indices.")
		indices = indices[indices >= 1 & indices < n]
	}
	
	if (length(indices) == 0)
	{
		return()
	}
	
	fit_summary = object@fits@summary
		
	feature_perm = 1:n
	if (sort_features)
	{
		feature_perm = c((1:n)[fit_summary$passed][order(-fit_summary$margin[fit_summary$passed])], (1:n)[!fit_summary$passed])
	}
	
	indices = feature_perm[indices]
	y = object@parameters@y

	for (i in indices)
	{
		grid.newpage()
		
		x = object@parameters@x[i,]
		feature = object@parameters@features[i]

		threshold = fit_summary$threshold[i]
		margin = fit_summary$margin[i]

		obj_plot = messinaSurvObjPlot(object, i) + ggtitle(sprintf("Objective Function: %s", object@parameters@perf_requirement$objective_type))
		
		if (bootstrap_type == "ci")			{ bootstrap_string = sprintf("Shaded area: %.0f%% CI", bootstrap_ci*100) }
		else if (bootstrap_type == "stdev")	{ bootstrap_string = "Shaded area: +/- 1 SD" }
		else 								{ bootstrap_string = "" }
		
		if (is.na(threshold))
		{
			km_plot_all = messinaSurvKMplotSingleGroup(y = y, bootstrap_type = bootstrap_type, bootstrap_ci = bootstrap_ci, nboot = nboot, parallel = parallel) + 
				ggtitle(paste("KM of Full Cohort (no valid threshold found)", bootstrap_string, sep = "\n"))

			pushViewport(viewport(layout = grid.layout(2, 1)))
			print(obj_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
			print(km_plot_all, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
			#popViewport(2)
		}
		else
		{
			km_plot_threshold = messinaSurvKMplot(y = y, group = (x > threshold)*1, bootstrap_type = bootstrap_type, bootstrap_ci = bootstrap_ci, nboot = nboot, parallel = parallel) + 
				ggtitle(paste("Separation at Threshold", bootstrap_string, sep = "\n")) + 
				theme(legend.position = "bottom")
			km_plot_lower_margin = messinaSurvKMplot(y = y, group = (x > threshold - margin/2)*1, bootstrap_type = bootstrap_type, bootstrap_ci = bootstrap_ci, nboot = nboot, parallel = parallel) + 
				ggtitle(paste("Separation at Lower Margin", bootstrap_string, sep = "\n")) + 
				theme(legend.position = "bottom")
			km_plot_upper_margin = messinaSurvKMplot(y = y, group = (x > threshold + margin/2)*1, bootstrap_type = bootstrap_type, bootstrap_ci = bootstrap_ci, nboot = nboot, parallel = parallel) + 
				ggtitle(paste("Separation at Upper Margin", bootstrap_string, sep = "\n")) + 
				theme(legend.position = "bottom")

			pushViewport(viewport(layout = grid.layout(2, 2)))
			print(obj_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
			print(km_plot_threshold, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
			print(km_plot_lower_margin, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
			print(km_plot_upper_margin, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
			#popViewport(4)
		}
	}
}


#' @import ggplot2
messinaSurvObjPlot = function(object, i)
{
	Threshold = Objective = NULL		# To shut up an R CMD check note for the later use of these in ggplot

	parameters = object@parameters
	objective_type = parameters@perf_requirement$objective_type
	objective_min = parameters@perf_requirement$min_objective

	threshold = object@fits@summary$threshold[i]
	margin = object@fits@summary$margin[i]

	objective_surfaces = object@fits@objective_surfaces[[i]]
	plot_data = data.frame(Objective = objective_surfaces$objective, Threshold = objective_surfaces$cutoff)

	cutoff_frac_points = quantile(parameters@x[i,], probs = c(parameters@minimum_group_fraction, 1 - parameters@minimum_group_fraction))
	cutoff_frac_ok = (parameters@x[i,] >= cutoff_frac_points[1]) && (parameters@x[i,] <= cutoff_frac_points[2])
	
	plot_data$Colour = c("darkgrey", "black")[cutoff_frac_ok + 1]
	
#	theplot = ggplot(data = plot_data, mapping = aes(x = Threshold, y = Objective, color = Colour)) + 
	theplot = ggplot(data = plot_data, mapping = aes(x = Threshold, y = Objective)) + 
		geom_line(alpha = 0.5) + 
		geom_point()

	if (zapsmall(parameters@minimum_group_fraction) != 0)
	{
		theplot = theplot + 
			geom_vline(xintercept = cutoff_frac_points, lty = "dotted", alpha = 0.5)
	}

	if (objective_type %in% c("tau", "reltau"))
	{
		theplot = theplot + 
			coord_cartesian(ylim = c(0, 1)) + 
			geom_hline(yintercept = c(objective_min, 1-objective_min), lty = "dotted") +
			geom_hline(yintercept = 0.5, lty = "solid", alpha = 0.5)
	}
	else if (objective_type == "coxcoef")
	{	
		y_limit = max(c(objective_min, abs(plot_data$Objective)[cutoff_frac_ok]))
		theplot = theplot + 
			coord_cartesian(ylim = c(-y_limit, y_limit)*1.3) + 
			geom_hline(yintercept = c(objective_min, -objective_min), lty = "dotted") + 
			geom_hline(yintercept = 0, lty = "solid", alpha = 0.5)
	}

	if (!is.na(threshold))
	{
		theplot = theplot + 
			geom_vline(xintercept = threshold, lty = "solid", alpha = 1)
			
		if (!is.na(margin) && margin != 0)
		{
			theplot = theplot + 
				geom_vline(xintercept = c(threshold - margin/2, threshold + margin/2), lty = "dashed", alpha = 0.5)
		}
	}
	
	return(invisible(theplot))
}


#' @importFrom survival survfit
#' @importFrom plyr llply
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


#' @importFrom plyr llply
calcBootstrapKaplanMeierEstimates = function(y, x, nboot, parallel)
{
	xvals = sort(as.character(unique(x)))
#	results1 = rlply(nboot, calcKaplanMeierEstimatesOnBootstrapSample(y, x))
	results1 = llply(as.list(1:nboot), function(dummy) calcKaplanMeierEstimatesOnBootstrapSample(y, x), .parallel = parallel)		# Using instead of rlply so we can get .parallel
	results2 = llply(xvals, function(xi) llply(results1, function(ri) ri[[xi]]))
	names(results2) = xvals
	return(results2)
}


#' @importFrom plyr llply laply
calcBootstrapKaplanMeierEstimatesAtTimes = function(y, x, nboot, times = sort(unique(y[,1])), parallel)
{
	ests = calcBootstrapKaplanMeierEstimates(y = y, x = x, nboot = nboot, parallel = parallel)
	ests_at_times = llply(ests, 
		function(ests_for_x) 
		{
			ests_at_time = laply(ests_for_x, 
				function(est) 
				{
					if (all(is.na(est$time)) || all(is.na(est$surv)))
					{
						return(rep(NA, length(times)))
					}
					else
					{
						return(approx(est$time, est$surv, times, method = "constant", f = 0, rule = 1, ties = "ordered")$y)
					}
				}
			)
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
	return(ests_at_times)
}


#' @import ggplot2
#' @importFrom plyr aaply
messinaSurvKMplotSingleGroup = function(y, bootstrap_type, bootstrap_ci, nboot, parallel)
{
	Time = Survival = NULL		# To shut up an R CMD check note for the later use of these in ggplot

	group = factor(rep("GROUP", nrow(y)))
	full_km = calcKaplanMeierEstimates(y, group)[[1]]
	
	if (bootstrap_type != "none")
	{
		boot_times = sort(unique(y[,1]))
		boot_km = calcBootstrapKaplanMeierEstimatesAtTimes(y = y, x = group, nboot = nboot, times = boot_times, parallel = parallel)[[1]]
	
		if (bootstrap_type == "ci")
		{
			boot_km_c = aaply(boot_km[,-1,drop=FALSE], 1, median, na.rm = TRUE)
			boot_km_l = aaply(boot_km[,-1,drop=FALSE], 1, quantile, probs = 1-bootstrap_ci, na.rm = TRUE)
			boot_km_u = aaply(boot_km[,-1,drop=FALSE], 1, quantile, probs = bootstrap_ci, na.rm = TRUE)
		}
		else if (bootstrap_type == "stdev")
		{
			boot_km_sd = aaply(boot_km[,-1,drop=FALSE], 1, sd, na.rm = TRUE)
			boot_km_c = rowMeans(boot_km[,-1])
			boot_km_l = boot_km_c - boot_km_sd
			boot_km_u = boot_km_c + boot_km_sd
		}
		else	{ stop("Unknown bootstrap_type") }
		
		boot_data = data.frame(	Time = boot_times, Survival = boot_km_c, SurvMin = boot_km_l, SurvMax = boot_km_u)
		boot_data = boot_data[!is.na(boot_data$Survival) & !is.na(boot_data$SurvMin) & !is.na(boot_data$SurvMax),]
	
		boot_poly_data = data.frame(	Time = rep(c(boot_times[-1], rev(boot_times)[-1]), each = 2),
										Survival = c(rep(boot_km_u[-1], each = 2)[-1], rep(rev(boot_km_l), each = 2)[-(length(boot_times)*2)]))
		boot_poly_data = boot_poly_data[!is.na(boot_poly_data$Survival),]
	}
	
	full_data = data.frame(	Time = full_km$time, Survival = full_km$surv)
	
	theplot = ggplot(data = full_data, aes(x = Time, y = Survival)) +
		geom_step(direction = "hv") + 
		coord_cartesian(ylim = c(0, 1))
	
	if (bootstrap_type != "none")
	{
		theplot = theplot + 
			geom_polygon(data = boot_poly_data, mapping = aes(x = Time, y = Survival), alpha = 0.2, linetype = 0) + 
			geom_step(data = boot_data, direction = "hv", alpha = 0.2)
	}
	
	return(invisible(theplot))
}


#' @import ggplot2
#' @importFrom plyr aaply llply
messinaSurvKMplot = function(y, group, bootstrap_type, bootstrap_ci, nboot, parallel)
{
	Time = Survival = Group = NULL		# To shut up an R CMD check note for the later use of these in ggplot

	stopifnot(all(group %in% c(TRUE, FALSE, 0, 1)))
	group = group*1
	
	full_km = calcKaplanMeierEstimates(y, group)
	
	if (bootstrap_type != "none")
	{
		boot_times = sort(unique(y[,1]))
		boot_km = calcBootstrapKaplanMeierEstimatesAtTimes(y = y, x = group, nboot = nboot, times = boot_times, parallel = parallel)
		names(boot_km) = c("0", "1")
	
		if (bootstrap_type == "ci")
		{
			boot_km_c = llply(boot_km, function(boots) aaply(boots[,-1,drop=FALSE], 1, median, na.rm = TRUE))
			boot_km_l = llply(boot_km, function(boots) aaply(boots[,-1,drop=FALSE], 1, quantile, probs = 1-bootstrap_ci, na.rm = TRUE))
			boot_km_u = llply(boot_km, function(boots) aaply(boots[,-1,drop=FALSE], 1, quantile, probs = bootstrap_ci, na.rm = TRUE))
		}
		else if (bootstrap_type == "stdev")
		{
			boot_km_sd = llply(boot_km, function(boots) aaply(boots[,-1,drop=FALSE], 1, sd, na.rm = TRUE))
			boot_km_c = llply(c("0", "1"), function(group) rowMeans(boot_km[[group]][,-1]))
			names(boot_km_c) = c("0", "1")
			boot_km_l = llply(c("0", "1"), function(group) boot_km_c[[group]] - boot_km_sd[[group]])
			boot_km_u = llply(c("0", "1"), function(group) boot_km_c[[group]] + boot_km_sd[[group]])
			names(boot_km_l) = c("0", "1")
			names(boot_km_u) = c("0", "1")
		}
		else	{ stop("Unknown bootstrap_type") }
		
		boot_data = data.frame(	Time = rep(boot_times, 2), 
								Survival = c(boot_km_c[["0"]], boot_km_c[["1"]]),
								SurvMin = c(boot_km_l[["0"]], boot_km_l[["1"]]),
								SurvMax = c(boot_km_u[["0"]], boot_km_u[["1"]]),
								Group = rep(c("<= Threshold", "> Threshold"), each = length(boot_times)))
		boot_data = boot_data[!is.na(boot_data$Survival) & !is.na(boot_data$SurvMin) & !is.na(boot_data$SurvMax),]
	
		boot_poly_data = data.frame(	Time = rep(rep(c(boot_times[-1], rev(boot_times)[-1]), each = 2), 2),
										Survival = c(c(rep(boot_km_u[["0"]][-1], each = 2)[-1], rep(rev(boot_km_l[["0"]]), each = 2)[-(length(boot_times)*2)]),
													 c(rep(boot_km_u[["1"]][-1], each = 2)[-1], rep(rev(boot_km_l[["1"]]), each = 2)[-(length(boot_times)*2)])),
										Group = rep(c("<= Threshold", "> Threshold"), each = (length(boot_times)-1)*4))
		boot_poly_data = boot_poly_data[!is.na(boot_poly_data$Survival),]
	}
	
	full_data = data.frame(	Time = c(full_km[["0"]]$time, full_km[["1"]]$time), 
							Survival = c(full_km[["0"]]$surv, full_km[["1"]]$surv), 
							Group = factor(c(rep("<= Threshold", nrow(full_km[["0"]])), rep("> Threshold", nrow(full_km[["1"]])))))
	
	theplot = ggplot(data = full_data, aes(x = Time, y = Survival, group = Group, col = Group)) +
		geom_step(direction = "hv") + 
		coord_cartesian(ylim = c(0, 1))
	
	if (bootstrap_type != "none")
	{
		theplot = theplot + 
			geom_polygon(data = boot_poly_data, mapping = aes(x = Time, y = Survival, group = Group, fill = Group), alpha = 0.2, linetype = 0) + 
			geom_step(data = boot_data, direction = "hv", alpha = 0.2)
	}
	
	return(invisible(theplot))
}
