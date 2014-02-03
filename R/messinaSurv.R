
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
		tau = (agree+tied/2)/(agree+disagree+tied)

		return(tau)
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


messinaSurvDoesObjectivePass = function(cutoff_obj, obj_min, obj_func)
{
	if (obj_func == "tau")
	{
		cutoff_obj_good_positive = (cutoff_obj >= obj_min) & (!is.na(cutoff_obj))
		cutoff_obj_good_negative = ((1-cutoff_obj) >= obj_min) & (!is.na(cutoff_obj))
	}
	else if (obj_func == "coxcoef")
	{
		cutoff_obj_good_positive = (cutoff_obj >= obj_min) & (!is.na(cutoff_obj))
		cutoff_obj_good_negative = ((-cutoff_obj) >= obj_min) & (!is.na(cutoff_obj))
	}
	else
	{
		stop(sprintf("Unknown objective function \"%s\"", obj_func))
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


messinaSurvTrain = function(x1, y, obj_min, obj_func = "tau")
{
	# Measure the concordance (resubstituted) on all cutoffs
	# possible in the training set.
	x_sorted = rle(sort(x1))$values
	x_cutoffs = (x_sorted[-1] + x_sorted[-length(x_sorted)])/2
	cutoff_obj = sapply(x_cutoffs, function(cutoff) messinaSurvObjectiveFunc((x1 > cutoff)*1, y, obj_func))
	
	# Find the widest region with a concordance index at least
	# as good as min_gamma.  Need to do this twice, to account for
	# CIs < 0.5, which could be good if the direction of the
	# classifier were reversed.
	cutoff_obj_good = messinaSurvDoesObjectivePass(cutoff_obj, obj_min, obj_func)
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
	
	result = list(threshold = best_threshold, margin = best_margin, direction = best_direction, cutoffs = x_cutoffs, obj = cutoff_obj, obj_func = obj_func, obj_min = obj_min)
	return(result)
}


messinaSurvSingleXDraw = function(x1, y, obj_min, obj_func, n_test)
{
	# TODO: Make sure there are enough events in the test set to be useful
	test_indices = sample.int(length(x1), n_test)
	train_indices = setdiff(1:length(x1), test_indices)
	x_test = x1[test_indices]
	x_train = x1[train_indices]
	y_test = y[test_indices,]
	y_train = y[train_indices,]

	# Train a MessinaSurv classifier on the training set
	fit = messinaSurvTrain(x_train, y_train, obj_min, obj_func)
	
	# Assess the classifier on the test set
	if (is.na(fit$threshold))
	{
		test_obj = NA
	}
	else
	{
		test_class = 1*(x_test*fit$direction > fit$threshold*fit$direction)
		test_obj = messinaSurvObjectiveFunc(test_class, y_test, obj_func)
	}
	return(test_obj)
}


messinaSurvSingleX = function(x1, y, obj_min, obj_func, n_boot, n_test)
{
	results = replicate(n_boot, messinaSurvSingleXDraw(x1, y, obj_min, obj_func, n_test))
	return(as.vector(results))
}


messinaSurvTrainOnSubset = function(x, y, obj_min, obj_func, subset, parallel)
{
	fits = llply(1:nrow(x), 
		function(xi) 
		{ 
			if (subset[xi])
			{
				fit = messinaSurvTrain(x[xi,], y, obj_min, obj_func)
				return(fit)
			}
			else
			{
				return(NULL)
			}
		}, .parallel = parallel, .progress = "time")
	
	return(fits)
}


#' Run the MessinaSurv algorithm to find features (eg. genes) that can define groups
#' of patients with very different survival times.
#'
#' TODO
#' 
#' @param x feature expression values, either supplied as an ExpressionSet, or as
#'   an object that can be converted to a matrix by as.matrix.  In the latter case,
#'   features should be in rows and samples in columns, with feature names taken
#'   from the rows of the object.
#' @param y a Surv object containing survival times and censoring status for each
#    sample in x.
#' @param obj_min the minimum acceptable sensitivity that a classifier separating
#'   the two groups of y must achieve.
#' @param obj_func the minimum acceptable specificity that a classifier separating
#'   the two groups of y must achieve.
#' @param f_train the fraction of samples to be used in the training splits of the
#'   bootstrap rounds.
#' @param n_boot the number of bootstrap rounds to use.
#' @param seed an optional random seed for the analysis.  If NULL, the R PRNG us used
#'   as-is.
#' @return an object of class "MessinaResult" containing the results of the analysis.
#' @export
#' @seealso \code{\link{MessinaResult-class}}
#' @seealso \code{\link{ExpressionSet}}
#' @references TODO
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
messinaSurv = function(x, y, obj_min, obj_func = "tau", f_train = 0.8, n_boot = 50, seed = NULL, parallel = FALSE)
{
	stopifnot(class(y) == "Surv")

	if (class(x) == "ExpressionSet")
	{
		features = featureNames(x)
		x = exprs(x)
	}
	else
	{
		x = as.matrix(x)
		features = rownames(x)
	}

	if (!is.null(seed))
	{
		set.seed(seed)
	}

	n_test = ceiling((1 - f_train) * nrow(y))

	cat("Performance bootstrapping...\n")
	boot_objs = aaply(x, 1, function(xi) messinaSurvSingleX(xi, y = y, obj_min = obj_min, obj_func = obj_func, n_boot = n_boot, n_test = n_test), .parallel = parallel, .progress = "time")

	rownames(boot_objs) = rownames(x)
	colnames(boot_objs) <- NULL
	
	# If fitting failed, set performance to the no-information value
	temp = boot_objs
	temp[is.na(temp)] = c("tau" = 0.5, "coxcoef" = 0)[obj_func]
	
	# Find the rows of x (if any) that passed the performance criterion on
	# the CV results.
	mean_obj = rowMeans(temp)
	obj_passes = mean_obj >= obj_min
	
	# Train MessinaSurv on just those rows of x that passed.
	cat("Final training...\n")
#	fits = messinaSurvTrainOnSubset(x, y, obj_min, obj_func, obj_passes, parallel = parallel)
	fits = messinaSurvTrainOnSubset(x, y, obj_min, obj_func, obj_passes | TRUE, parallel = parallel)
	
	thresholds = sapply(fits, function(f) ifelse(is.null(f), NA, f$threshold))
	posks = sapply(fits, function(f) ifelse(is.null(f), NA, f$direction)) == 1
	margins = sapply(fits, function(f) ifelse(is.null(f), NA, f$margin))
	
    result = list(  problem = "survival",
                    parameters = list(  x = x, 
                                        y = y, 
                                        features = features, 
                                        perf_requirement = list(obj_func = obj_func, obj_min = obj_min),
                                        f_train = f_train, 
                                        n_boot = n_boot, 
                                        seed = seed),
                    classifier = list(  type = as.factor(ifelse(obj_passes, "Threshold", NA)), 
                                        threshold = thresholds, 
                                        ptrue = rep(NA, nrow(x)),
                                        posk = posks), 
					fits = fits,
                    margin = margins, 
                    psuccessful = rowMeans(!is.na(boot_objs)),
                    passed = obj_passes,
                    perf = list(mean = mean_obj, var = NA),
                    bootstraps = boot_objs)
    class(result) = "MessinaResult"
	
	return(result)
}


