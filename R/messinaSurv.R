# messinaSurv.R: Main Messina survival fitting functions
# 
# Copyright 2014 Mark Pinese
#
# This file is distributed under the terms of the Eclipse Public 
# License v1.0, available at:
# https://www.eclipse.org/org/documents/epl-v10.html


#' Find optimal prognostic features using the Messina algorithm
#' 
#' Run the MessinaSurv algorithm to find features (eg. genes) that can define groups
#' of patients with very different survival times.
#' 
#' The MessinaSurv algorithm aims to identify features for which patients with high signal
#' and patients with low signal have very different survival outcomes.  This is achieved by 
#' definining an objective function which assigns a numerical value to how strongly the 
#' survival in two groups of patients differs, then assessing the value of this objective
#' at different signal levels of each feature.  Those features for which, at a given signal
#' level, the objective function is consistently above a user-supplied minimum level, are
#' selected by MessinaSurv as being single-feature survival predictors.
#' 
#' MessinaSurv has applications as an algorithm to identify features that are survival-related,
#' as well as a principled method to identify threshold signal values to separate a cohort into
#' poor- and good-prognosis subgroups.  It can also be used as a feature filter, selecting and
#' discretising survival-related features before they are input into a multivariate predictor.
#' 
#' @section Objective functions: MessinaSurv uses the value of its objective function as a
#'   measure of the strength of the difference in survival of the two patient groups 
#'   defined by the threshold.  Three objective functions are currently defined:
#'   \describe{
#'     \item{"coxcoef"}{The coefficient of a Cox proportional hazards fit to the model Surv ~ I(x > T),
#'     where x is the feature signal level, and T is the threshold being tested.  Range is (-inf, inf),
#'     with a no-information value of 0;
#'     positive values indicate that the subgroup defined by signal above the threshold fails sooner.}
#'     \item{"tau"}{Kendall's tau for survival data, defined as (concordant + tied/2) / (concordant + discordant + tied),
#'       where concordant is the number of concordant group/survival pairs, discordant is
#'       the number of discordant group/survival pairs, and tied is the total number of tied pairs,
#'       counting both group and survival ties.  Concordance is calculated expecting that
#'       samples with signal exceeding the threshold will fail sooner.  Range is [0, 1], with
#'       a no-information value of 0.5.  Note that the ties terms naturally penalize very high or low 
#'       thresholds, and so this objective is inappropriate if somewhat unbalanced subgroups are expected 
#'       to be present in the data.}
#'     \item{"reltau"}{tau, normalized to remove the ties penalty.  Defined as agree / (agree + disagree).
#'       Range is [0, 1], with a no-information value of 0.5.  Although the ties penalty of tau
#'       is removed, and this method is thus suitable for finding unbalanced subgroups, it is
#'       now unstable at extreme threshold values (as in these cases, agree + disagree -> 0).
#'       For this reason, min_group_frac must be set to a modest value when using "reltau", to preserve stability. }
#'   }
#'   Methods "coxcoef" and "reltau" show instability for very high and low threshold values,
#'   and so should be used with an appropriate value of min_group_frac for stable fits.  Method
#'   "tau" is stable to extreme threshold values, and therefore will tolerate min_group_frac = 0,
#'   however note that "tau" naturally penalizes small subgroups, and is therefore a poor
#'   choice unless you wish to find approximately equal-sized subgroups.
#' 
#' @section Minimum group fraction: The parameter min_group_frac limits the size of
#'   the smallest subgroups that messinaSurv can select.  As the groups become smaller,
#'   the "reltau" and "coxcoef" objective functions become unstable, and can generate
#'   spurious results.  These are seen on the diagnostics produced by the messina
#'   plot functions as very high objective values at very low and high threshold values.
#'   To control these results, set min_group_frac to a high enough value that the
#'   objective functions reliably fit.  Generally, max(0.1, 10/N), where N is the total
#'   number of patients, is sufficient.  Keep in mind that setting this parameter too high 
#'   will limit messinaSurv's ability to identify small subsets of patients with dramatically
#'   different survival from the rest: the smallest subset that will be reliably identified
#'   is min_group_frac of patients.
#' 
#' @param x feature expression values, either supplied as an ExpressionSet, or as
#'   an object that can be converted to a matrix by as.matrix.  In the latter case,
#'   features should be in rows and samples in columns, with feature names taken
#'   from the rows of the object.
#' @param y a Surv object containing survival times and censoring status for each
#    sample in x.
#' @param obj_min the minimum acceptable value of the objective metric.  The metric used
#'   is specified by the parameter obj_func.
#' @param obj_func the metric function that measures the difference in survival between
#'   patients with feature values above, and below, the threshold.  Valid values are 
#'   "tau", "reltau", or "coxcoef"; see details for more information.
#' @param min_group_frac the size of the smallest sample group that is allowed to be
#'   generated by thresholding, as a fraction of the total sample.  The default value
#'   of 0.1 means that no thresholds will be selected that result in a sample split
#'   yielding a group of smaller than 10% of the samples.  A modest value of this
#'   parameter increases the stability of the "reltau" and "coxcoef" objectives, which
#'   tend to become unstable as the number of samples in a group becomes very low; see
#'   details.
#' @param f_train the fraction of samples to be used in the training splits of the
#'   bootstrap rounds.
#' @param n_boot the number of bootstrap rounds to use.
#' @param seed an optional random seed for the analysis.  If NULL, the R PRNG is used
#'   as-is.
#' @param parallel should calculations be parallelized using the doMC framework?  If NULL, 
#'   parallel mode is used if the doMC library is loaded, and more than one
#'   core has been registered with registerDoMC().  Note that no progress bar is
#'   displayed in parallel mode.
#' @param silent be completely silent (except for error and warning messages)?
#' 
#' @return an object of class "MessinaSurvResult" containing the results of the analysis.
#' 
#' @importFrom plyr aaply llply
#'
#' @export
#' @seealso \code{\link{MessinaSurvResult-class}}
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @seealso \code{\link{messina}}
#' @seealso \code{\link{messinaDE}}
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
#' 
#' @examples
#' ## Load a subset of the TCGA renal clear cell carcinoma data
#' ## as an example.
#' data(tcga_kirc_example)
#' 
#' ## Run the messinaSurv analysis on these data.  Use a tau
#' ## objective, with a minimum performance of 0.6.  Note that
#' ## messinaSurv analyses are very computationally-intensive,
#' ## so in actual use multicore use with doMC and parallel = TRUE
#' ## is strongly recommended.
#' fit = messinaSurv(kirc.exprs, kirc.surv, obj_func = "tau", obj_min = 0.6)
#'
#' fit
#' plot(fit)
messinaSurv = function(x, y, obj_min, obj_func, min_group_frac = 0.1, f_train = 0.8, n_boot = 50, seed = NULL, parallel = NULL, silent = FALSE)
{
	stopifnot(class(y) == "Surv")
	
	if (is.null(parallel))
	{
		if (("doMC" %in% .packages()) && require(foreach))
		{
			parallel = foreach::getDoParWorkers() > 1
		}
		else
		{
			parallel = FALSE
		}
	}

	if (class(x) == "ExpressionSet")
	{
		if (require(Biobase) == FALSE)
		{
			stop("Bioconductor package Biobase must be available to use data in ExpressionSet objects")
		}
		features = Biobase::featureNames(x)
		samples = Biobase::sampleNames(x)
		x = Biobase::exprs(x)
	}
	else
	{
		x = as.matrix(x)
		features = rownames(x)
		samples = colnames(x)
		
		if (is.null(features))	{ features = paste("F", 1:nrow(x), sep = "") }
		if (is.null(samples))	{ samples = paste("S", 1:ncol(x), sep = "") }
	}
	
	n_boot = as.integer(n_boot)
	seed = as.integer(seed)

	if (length(seed) != 0)
	{
		set.seed(seed)
	}

	n_test = ceiling((1 - f_train) * nrow(y))

	if (!silent) { message("Performance bootstrapping...") }
	boot_objs = aaply(x, 1, function(xi) messinaSurvSingleX(xi, y = y, min_group_frac, obj_min = obj_min, obj_func = obj_func, n_boot = n_boot, n_test = n_test), .parallel = parallel, .progress = ifelse(parallel || silent, "none", "time"))

	rownames(boot_objs) = rownames(x)
	colnames(boot_objs) <- NULL
	
	# If fitting failed, set performance to the no-information value
	temp = boot_objs
	temp[is.na(temp)] = c("tau" = 0.5, "reltau" = 0.5, "coxcoef" = 0)[obj_func]
	
	# Find the rows of x (if any) that passed the performance criterion on
	# the CV results.
	mean_obj = rowMeans(temp)
	obj_passes = mean_obj >= obj_min
	
	# Train MessinaSurv on all rows of x (commented original was only on 
	# passing rows, but it's proven useful to have all fit trajectories
	# for diagnostics, and it adds little extra compute time above and
	# beyond the LOOCV).
	if (!silent) { message("Final training...") }
#	fits = messinaSurvTrainOnSubset(x, y, min_group_frac, obj_min, obj_func, obj_passes, parallel = parallel, silent = silent)
	fits = messinaSurvTrainOnSubset(x, y, min_group_frac, obj_min, obj_func, obj_passes | TRUE, parallel = parallel, silent = silent)
	
	thresholds = sapply(fits, function(f) ifelse(is.null(f), NA, f$threshold))
	posks = sapply(fits, function(f) ifelse(is.null(f), NA, f$direction)) == 1
	margins = sapply(fits, function(f) ifelse(is.null(f), NA, f$margin))
	
	parameters = .MessinaParameters(x = x,
									y = y,
									features = features,
									samples = samples,
									perf_requirement = list(objective_type = obj_func,
															min_objective = obj_min),
									minimum_group_fraction = min_group_frac,
									training_fraction = f_train,
									num_bootstraps = n_boot,
									prng_seed = seed)
	
	objective_surfaces = llply(fits, function(fit) data.frame(cutoff = fit$cutoffs, objective = fit$obj))
	
	fits2 = .MessinaFits(summary = data.frame(	passed = obj_passes,
												type = as.factor(ifelse(obj_passes, "Threshold", NA)),
												threshold = thresholds,
												posk = posks,
												margin = margins,
												ptrue = rep(NA, nrow(x)),
												psuccessful = rowMeans(!is.na(boot_objs))),
						objective_surfaces = objective_surfaces)
	
	result = .MessinaSurvResult(problem_type = "survival",
								parameters = parameters,
								perf_estimates = data.frame(mean_obj = mean_obj),
								fits = fits2)
	
	return(result)
}


#' @importFrom survival survConcordance.fit coxph
messinaSurvObjectiveFunc = function(x, y, func)
{
	if (func == "tau")
	{
		counts = survival::survConcordance.fit(y, x)
		agree = counts["concordant"]
		disagree = counts["discordant"]
		tied.time = counts["tied.time"]
		tied.risk = counts["tied.risk"]
		tied = tied.time + tied.risk
		tau = (agree+tied/2)/(agree+disagree+tied)

		return(tau)
	}
	else if (func == "reltau")
	{
		counts = survival::survConcordance.fit(y, x)
		agree = counts["concordant"]
		disagree = counts["discordant"]
		reltau = agree/(agree+disagree)

		return(reltau)
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
	if (length(cutoff_obj) == 0)
	{
		return(list(pos = logical(0), neg = logical(0)))
	}

	if (obj_func == "tau" || obj_func == "reltau")
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


messinaSurvTrain = function(x1, y, min_group_frac, obj_min, obj_func)
{
	# Measure the concordance (resubstituted) on all cutoffs
	# possible in the training set.
	x_sorted = rle(sort(x1))$values
	if (length(x_sorted) == 1)	{ x_cutoffs = x_sorted } else { x_cutoffs = (x_sorted[-1] + x_sorted[-length(x_sorted)])/2 }
	x_cutoffs = c(-Inf, x_cutoffs, Inf)
	x_cutoffs = x_cutoffs[!duplicated(x_cutoffs)]
	cutoff_obj = sapply(x_cutoffs, function(cutoff) messinaSurvObjectiveFunc((x1 > cutoff)*1, y, obj_func))
	
	# Find the widest region with a concordance index at least
	# as good as min_gamma.  Need to do this twice, to account for
	# CIs < 0.5, which could be good if the direction of the
	# classifier were reversed.
	cutoff_frac = sapply(x_cutoffs, function(cutoff) mean(x1 > cutoff))
	cutoff_frac = pmin(cutoff_frac, 1 - cutoff_frac)
	cutoff_frac_good = cutoff_frac >= min_group_frac

	cutoff_obj_good = messinaSurvDoesObjectivePass(cutoff_obj, obj_min, obj_func)
	cutoff_obj_good_positive = cutoff_obj_good$pos & cutoff_frac_good
	cutoff_obj_good_negative = cutoff_obj_good$neg & cutoff_frac_good
	
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
	
	result = list(threshold = best_threshold, margin = best_margin, direction = best_direction, cutoffs = x_cutoffs, obj = cutoff_obj, min_group_frac = min_group_frac, obj_func = obj_func, obj_min = obj_min)
	return(result)
}


messinaSurvSingleXDraw = function(x1, y, min_group_frac, obj_min, obj_func, n_test)
{
	# TODO: Make sure there are enough events in the test set to be useful
	test_indices = sample.int(length(x1), n_test)
	train_indices = setdiff(1:length(x1), test_indices)
	x_test = x1[test_indices]
	x_train = x1[train_indices]
	y_test = y[test_indices,]
	y_train = y[train_indices,]

	# Train a MessinaSurv classifier on the training set
	fit = messinaSurvTrain(x_train, y_train, min_group_frac, obj_min, obj_func)
	
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


messinaSurvSingleX = function(x1, y, min_group_frac, obj_min, obj_func, n_boot, n_test)
{
	results = replicate(n_boot, messinaSurvSingleXDraw(x1, y, min_group_frac, obj_min, obj_func, n_test))
	return(as.vector(results))
}


#' @importFrom plyr aaply llply
messinaSurvTrainOnSubset = function(x, y, min_group_frac, obj_min, obj_func, subset, parallel, silent)
{
	fits = llply(1:nrow(x), 
		function(xi) 
		{ 
			if (subset[xi])
			{
				fit = messinaSurvTrain(x[xi,], y, min_group_frac, obj_min, obj_func)
				return(fit)
			}
			else
			{
				return(NULL)
			}
		}, .parallel = parallel, .progress = ifelse(parallel || silent, "none", "time"))
	
	return(fits)
}

