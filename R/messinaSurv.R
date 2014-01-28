
messinaSurvObjectiveFunc = function(x, y, func)
{
	if (func == "tau")
	{
		counts = survival:::survConcordance.fit(y, x)
		agree = counts["concordant"]
		disagree = counts["discordant"]
		tied.time = counts["tied.time"]
		tied.risk = counts["tied.risk"]
		tied = tied.time + tied.risk
		return((agree+tied/2)/(agree+disagree+tied))
	}
	else if (func == "coxcoef")
	{
		fit = try(coxph(y ~ x))
		if (typeof(fit) == "try-error")	{ return(NA) }
		return(coef(fit))
	}
	else
	{
		stop(sprintf("Unknown objective function \"%s\"", func))
	}
}


messinaSurvDoesObjectivePass = function(cutoff_obj, obj.min, obj.func)
{
	if (obj.func == "tau")
	{
		cutoff_obj_good_positive = (cutoff_obj >= obj.min) & (!is.na(cutoff_obj))
		cutoff_obj_good_negative = ((1-cutoff_obj) >= obj.min) & (!is.na(cutoff_obj))
	}
	else if (obj.func == "coxcoef")
	{
		cutoff_obj_good_positive = (cutoff_obj >= obj.min) & (!is.na(cutoff_obj))
		cutoff_obj_good_negative = ((-cutoff_obj) >= obj.min) & (!is.na(cutoff_obj))
	}
	else
	{
		stop(sprintf("Unknown objective function \"%s\"", obj.func))
	}
	
	return(list(pos = cutoff_obj_good_positive, neg = cutoff_obj_good_negative))
}


messinaSurvFindBestThreshold = function(x, f)
{
	if (sum(f) == 0)
	{
		result = list(margin = 0, threshold = NA)
	}
	else
	{
		region_starts = which(c(f[1], f[-1] & (!f[-length(f)])))
		region_ends = which(c(f[-length(f)] & (!f[-1]), f[length(f)]))
		region_margins = x[region_ends] - x[region_starts]
		region_thresholds = x[region_starts] + region_margins/2
		result = list(margin = max(region_margins), threshold = region_thresholds[which.max(region_margins)])
	}
	return(result)
}


messinaSurvTrain = function(x1, y, obj.min, obj.func = "tau")
{
	# Measure the concordance (resubstituted) on all cutoffs
	# possible in the training set.
	x_sorted = rle(sort(x1))$values
	x_cutoffs = (x_sorted[-1] + x_sorted[-length(x_sorted)])/2
	cutoff_obj = sapply(x_cutoffs, function(cutoff) messinaSurvObjectiveFunc((x1 > cutoff)*1, y, obj.func))
	
	# Find the widest region with a concordance index at least
	# as good as min_gamma.  Need to do this twice, to account for
	# CIs < 0.5, which could be good if the direction of the
	# classifier were reversed.
	cutoff_obj_good = messinaSurvDoesObjectivePass(cutoff_obj, obj.min, obj.func)
	cutoff_obj_good_positive = cutoff_obj_good$pos
	cutoff_obj_good_negative = cutoff_obj_good$neg
	
	# Find the best threshold and direction
	best_pos = messinaSurvFindBestThreshold(x_sorted, cutoff_obj_good_positive)
	best_neg = messinaSurvFindBestThreshold(x_sorted, cutoff_obj_good_negative)
	if (best_pos$margin > best_neg$margin)
	{
		best_direction = 1
		best_threshold = best_pos$threshold
		best_margin = best_pos$margin
	}
	else
	{
		best_direction = -1
		best_threshold = best_neg$threshold
		best_margin = best_neg$margin
	}
	
	if (is.na(best_threshold))
	{
		best_direction = NA
	}
	
	result = list(threshold = best_threshold, margin = best_margin, direction = best_direction, cutoffs = x_cutoffs, obj = cutoff_obj, obj.func = obj.func, obj.min = obj.min)
	return(result)
}


messinaSurvSingleXDraw = function(x1, y, obj.min, obj.func, n_test)
{
	test_indices = sample.int(length(x1), n_test)
	train_indices = setdiff(1:length(x1), test_indices)
	x_test = x1[test_indices]
	x_train = x1[train_indices]
	y_test = y[test_indices,]
	y_train = y[train_indices,]

	# Train a MessinaSurv classifier on the training set
	fit = messinaSurvTrain(x_train, y_train, obj.min, obj.func)
	
	# Assess the classifier on the test set
	if (is.na(fit$threshold))
	{
		test_obj = NA
	}
	else
	{
		test_class = 1*(x_test*fit$direction > fit$threshold*fit$direction)
		test_obj = messinaSurvObjectiveFunc(test_class, y_test, obj.func)
	}
	return(test_obj)
}


messinaSurvSingleX = function(x1, y, obj.min, obj.func, n_boot, n_test)
{
	results = replicate(n_boot, messinaSurvSingleDraw(x1, y, obj.min, obj.func, n_test))
	return(results)
}


messinaSurvTrainOnSubset = function(x, y, obj.min, obj.func, subset, parallel)
{
	fits = adply(x, 
		function(xi) 
		{ 
			if (subset[xi])
			{
				fit = messinaSurvTrain(x[xi,], y, obj.min, obj.func)
				return(c(fit$threshold, fit$margin, fit$direction))
			}
			else
			{
				return(c(NA, NA, NA))
			}
		}, .parallel = parallel, .progress = "time")
	
	colnames(fits) = c("Threshold", "Margin", "Direction")
	return(fits)
}


messinaSurv = function(x, y, obj.min, obj.func = "tau", n_boot = 50, test_fraction = 0.2, seed = NULL, parallel = FALSE)
{
	stopifnot(class(y) == "Surv")
	
	if (!is.null(seed))
	{
		set.seed(seed)
	}

	n_test = ceiling(test_fraction * nrow(y))

	boot_objs = aaply(x, 1, function(xi) messinaSurvSingleX(xi, y = y, obj.min = obj.min, obj.func = obj.func, n_boot = n_boot, n_test = n_test), .parallel = parallel, .progress = "time")

	rownames(boot_objs) = rownames(x)
	colnames(boot_objs) <- NULL
	
	# If fitting failed, set performance to the no-information value
	temp = boot_objs
	temp[is.na(temp)] = c("tau" = 0.5, "coxcoef" = 0)[obj.func]
	
	# Find the rows of x (if any) that passed the performance criterion on
	# the CV results.
	mean_obj = rowMeans(temp)
	obj_passes = mean_obj >= obj.min
	
	# Train MessinaSurv on just those rows of x that passed.
	fits = messinaSurvTrainOnSubset(x, y, obj.min, obj.func, obj_passes)
	
	result = list(summary = data.frame(PassedObj = obj_passes, MeanObj = mean_obj, Threshold = fits$Threshold, Margin = fits$Margin, Direction = fits$Direction), boot.objs = boot_objs, settings = list(obj.min = obj.min, obj.func = obj.func, n_boot = n_boot, test_fraction = test_fraction, seed = seed))
	return(result)
}
