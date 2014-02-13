setMethod("plot", signature = signature(x = "MessinaClassResult", y = "missing"), definition = function(x, y, ...) messinaClassPlot(object = x, ...))

setMethod("plot", signature = signature(x = "MessinaSurvResult", y = "missing"), definition = function(x, y, ...) messinaSurvPlot(object = x, ...))


messinaClassPlot = function(object, indices = c(1), sort_features = TRUE, bootstrap_type = "none", bootstrap_ci = 0.90, nboot = ifelse(bootstrap_type == "ci", 50/(1-bootstrap_ci), 50))
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

	n = nrow(object@parameters@x)
	
	if (any(indices > n) || any(indices < 1))
	{
		warning("Warning: Some feature indices were out of range.  Dropping invalid indices.")
		indices = indices[indices >= 1 & indices < n]
	}
	
	if (length(indices) == 0)
	{
		return(invisible(NULL))
	}
	
	fit_summary = object@fits@summary
		
	feature_perm = 1:n
	if (sort_features)
	{
		feature_perm = c((1:n)[fit_summary$PassedPerf][order(-fit_summary$Margin[fit_summary$PassedPerf])], (1:n)[!fit_summary$PassedPerf])
	}
	
	indices = feature_perm[indices]
	y = object@parameters@y

	for (i in indices)
	{
		x = object@parameters@x[i,]
		feature = object@parameters@features[i]

		threshold = fit_summary$threshold[i]
		margin = fit_summary$margin[i]

		obj_plot = messinaSurvObjPlot(object, i) + ggtitle("Objective Function")
		
		if (bootstrap_type == "ci")			{ bootstrap_string = sprintf("Shaded area: %.0f%% CI", bootstrap_ci*100) }
		else if (bootstrap_type == "stdev")	{ bootstrap_string = "Shaded area: +/- 1 SD" }
		else 								{ bootstrap_string = "" }
		
		# TODO:
		if (is.na(threshold))
		{
			km_plot_all = messinaSurvKMplotSingleGroup(y, bootstrap_type, bootstrap_ci, nboot) + ggtitle(paste("KM of Full Cohort (no valid threshold found)", bootstrap_string, sep = "\n"))
			
			theplot = arrangeGrob(obj_plot, km_plot_all, main = sprintf("MessinaSurv Fit: Feature %s", feature))
		}
		else
		{
			km_plot_threshold = messinaSurvKMplot(y, (x > threshold)*1, bootstrap_type, bootstrap_ci, nboot) + ggtitle(paste("Separation at Threshold", bootstrap_string, sep = "\n")) + theme(legend.position = "bottom")
			km_plot_lower_margin = messinaSurvKMplot(y, (x > threshold - margin/2)*1, bootstrap_type, bootstrap_ci, nboot) + ggtitle(paste("Separation at Lower Margin", bootstrap_string, sep = "\n")) + theme(legend.position = "bottom")
			km_plot_upper_margin = messinaSurvKMplot(y, (x > threshold + margin/2)*1, bootstrap_type, bootstrap_ci, nboot) + ggtitle(paste("Separation at Upper Margin", bootstrap_string, sep = "\n")) + theme(legend.position = "bottom")
			
			theplot = arrangeGrob(obj_plot, km_plot_threshold, km_plot_lower_margin, km_plot_upper_margin, main = sprintf("MessinaSurv Fit: Feature %s", feature))
		}
		
		print(theplot)
	}
	
	invisible(NULL)
}


messinaSurvPlot = function(object, indices = c(1), sort_features = TRUE, plot_type = "point")
{
	#Sample = Value = Class = NULL		# To shut up an R CMD check note for the later use of these in ggplot
	
	n = nrow(object@parameters@x)
	
	if (any(indices > n) || any(indices < 1))
	{
		warning("Warning: Some feature indices were out of range.  Dropping invalid indices.")
		indices = indices[indices >= 1 & indices < n]
	}
	
	if (length(indices) == 0)
	{
		return(invisible(NULL))
	}
	
	fit_summary = object@fits@summary
	
	feature_perm = 1:n
	if (sort_features)
	{
		feature_perm = c((1:n)[fit_summary$PassedPerf][order(-fit_summary$Margin[fit_summary$PassedPerf])], (1:n)[!fit_summary$PassedPerf])
	}
	
	indices = feature_perm[indices]
	
	for (i in indices)
	{
		feature = object@parameters@features[i]
		x = result$parameters$x[i,]

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

		if (type == "point")	{ theplot = theplot + geom_point(stat = "identity") }
		else if (type == "bar")	{ theplot = theplot + geom_bar(stat = "identity", width = 0.5) }
		else 					{ stop(sprintf("Invalid Messina plot type: \"%s\"", type)) }	

		if (!is.na(threshold))
		{
			theplot = theplot + 
				geom_hline(yintercept = c(threshold - margin/2, threshold, threshold + margin/2), linetype = c("dashed", "solid", "dashed"), lwd = 0.7)
		}
		
		print(theplot)
	}
	
	invisible(NULL)
}



messinaSurvObjPlot = function(object, i)
{
	parameters = object@parameters
	objective_type = parameters@perf_requirement$objective_type
	objective_min = parameters@perf_requirement$min_objective

	threshold = object@fit@summary$threshold[i]
	margin = object@fit@summary$margin[i]

	objective_surface = object@fits@objective_surface[[i]]
	plot_data = data.frame(Objective = objective_surface$objective, Cutoff = objective_surface$cutoff)
	
	theplot = ggplot(data = plot_data, mapping = aes(x = Cutoff, y = Objective)) + 
		geom_line(alpha = 0.5) + 
		geom_point()
		
	if (objective_type == "tau")
	{
		theplot = theplot + 
			coord_cartesian(ylim = c(0, 1)) + 
			geom_hline(yintercept = c(objective_min, 1-objective_min), lty = "dotted") +
			geom_hline(yintercept = 0.5, lty = "solid", alpha = 0.5)
	}
	else if (objective_type == "coxcoef")
	{	
		theplot = theplot + 
			geom_hline(yintercept = c(objective_min, -objective_min), lty = "dotted") + 
			geom_hline(yintercept = 0, lty = "solid", alpha = 0.5)
	}

	if (!is.na(threshold))
	{
		theplot = theplot + 
			geom_vline(xintercept = threshold, lty = "solid", alpha = 1)
			
		if (margin != 0)
		{
			theplot = theplot + 
				geom_vline(xintercept = c(threshold - margin/2, threshold + margin/2), lty = "dashed", alpha = 0.5)
		}
	}
	
	theplot
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


messinaSurvKMplotSingleGroup = function(y, bootstrap_type, bootstrap_ci, nboot)
{
	group = factor(rep("GROUP", nrow(y)))
	full_km = calcKaplanMeierEstimates(y, group)[[1]]

	if (bootstrap_type != "none")
	{
		boot_times = sort(unique(y[,1]))
		boot_km = calcBootstrapKaplanMeierEstimatesAtTimes(y, group, nboot, boot_times)[[1]]

		if (bootstrap_type == "ci")
		{
			boot_km_c = aaply(boot_km[,-1,drop=FALSE], 1, median, na.rm = TRUE)
			boot_km_l = aaply(boot_km[,-1,drop=FALSE], 1, quantile, probs = 1-bootstrap_ci, na.rm = TRUE)
			boot_km_u = aaply(boot_km[,-1,drop=FALSE], 1, quantile, probs = bootstrap_ci, na.rm = TRUE)
		}
		else
		{
			boot_km_sd = aaply(boot_km[,-1,drop=FALSE], 1, sd, na.rm = TRUE)
			boot_km_c = rowMeans(boot_km[,-1])
			boot_km_l = boot_km_c - boot_km_sd
			boot_km_u = boot_km_c + boot_km_sd
		}
		
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

	theplot	
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
	
	theplot	
}
