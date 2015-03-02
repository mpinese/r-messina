# messinaSurv.R: Main Messina survival fitting functions
# 
# Copyright 2014-2015 Mark Pinese
#
# This file is distributed under the terms of the Eclipse Public 
# License v1.0, available at:
# https://www.eclipse.org/org/documents/epl-v10.html
#
# Changelog:
# 20150301	Replaced the original coxcoef implementation with a much faster and
#			more stable approximate one.  The approximation is very close, so
#			only erroneous results due to unstable fits should be affected.


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
#'     \item{"coxcoef"}{The approximate coefficient of a Cox proportional hazards fit to the model Surv ~ I(x > T),
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
#' @import foreach
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
messinaSurv = function(x, y, objective, min_group_frac = 0.1, f_train = 0.8, n_boot = 50, seed = NULL, parallel = NULL, silent = FALSE)
{
	stopifnot(class(y) == "Surv")
	
	if (is.null(parallel))
	{
		if ("doMC" %in% .packages())
		{
			parallel = getDoParWorkers() > 1
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
	boot_objs = llply(1:nrow(x), function(xi) 
		messinaSurvSingleX(x[xi,], y = y, min_group_frac, objective = objective, n_boot = n_boot, n_test = n_test), .parallel = parallel, .progress = ifelse(parallel || silent, "none", "time"))
	boot_objs_passed = laply(boot_objs, function(boot_obj) as.logical(boot_obj["pass_good_positive",]))
	boot_objs_metrics = laply(boot_objs, function(boot_obj) as.numeric(boot_obj["metric",]))

	median_obj_metric = apply(boot_objs_metrics, 1, median, na.rm = TRUE)

	rownames(boot_objs_passed) = rownames(x)
	colnames(boot_objs_passed) <- NULL
	
	# If fitting failed, set performance to the no-information value
	temp = boot_objs_passed
	temp[is.na(temp)] = FALSE
	
	# Find the rows of x (if any) that passed the performance criterion on
	# the CV results.
	obj_passes = rowMeans(temp) >= 0.5
	
	# Train MessinaSurv on all rows of x
	if (!silent) { message("Final training...") }
	fits = messinaSurvTrainOnSubset(x, y, min_group_frac, objective, obj_passes | TRUE, parallel = parallel, silent = silent)
	
	thresholds = sapply(fits, function(f) ifelse(is.null(f), NA, f$threshold))
	posks = sapply(fits, function(f) ifelse(is.null(f), NA, f$direction)) == 1
	margins = sapply(fits, function(f) ifelse(is.null(f), NA, f$margin))
	
	parameters = .MessinaParameters(x = x,
									y = y,
									features = features,
									samples = samples,
									perf_requirement = list(objective = objective),
									minimum_group_fraction = min_group_frac,
									training_fraction = f_train,
									num_bootstraps = n_boot,
									prng_seed = seed)
	
	objective_surfaces = llply(fits, function(fit) data.frame(cutoff = fit$cutoffs, objective = fit$objective_value))
	
	fits2 = .MessinaFits(summary = data.frame(	passed = obj_passes,
												type = as.factor(ifelse(obj_passes, "Threshold", NA)),
												threshold = thresholds,
												posk = posks,
												margin = margins,
												ptrue = rep(NA, nrow(x)),
												psuccessful = rowMeans(!is.na(boot_objs_metrics))),
						objective_surfaces = objective_surfaces)
	
	result = .MessinaSurvResult(problem_type = "survival",
								parameters = parameters,
								perf_estimates = data.frame(median_metric = median_obj_metric),
								fits = fits2)
	
	return(result)
}


#' @export
messinaSurvObj.CoxCoef = function(coxcoef_threshold)
{
	stopifnot(coxcoef_threshold > 0)

	f <- function(x, y)
	{
		sort_perm = order(y[,1])
		time = y[sort_perm,1]
		event = y[sort_perm,2]
		x = x[sort_perm]
		zstat = messinaSurvLRT(as.logical(x), time, event)
		coef = zstat * sqrt(4/sum(event))

		return(list(pass_good_positive = -coef >= coxcoef_threshold, pass_good_negative = coef >= coxcoef_threshold, metric = coef))
	}
	attr(f, "ObjCall") <- sprintf("messinaSurvObj.CoxCoef(coxcoef_threshold = %s)", as.character(coxcoef_threshold))
	attr(f, "PlotThresh") <- c(coxcoef_threshold, 0, -coxcoef_threshold)
	f
}


#' @export
messinaSurvObj.Tau = function(tau_threshold)
{
	stopifnot(tau_threshold > 0.5 && tau_threshold < 1)

	f <- function(x, y)
	{
		sort_perm = order(y[,1])
		time = y[sort_perm,1]
		event = y[sort_perm,2]
		x = x[sort_perm]
		counts = messinaSurvConcordance(as.logical(x), time, event)

		tau = (counts$concordant + counts$ties/2) / (counts$concordant + counts$discordant + counts$ties)

		return(list(pass_good_positive = tau >= tau_threshold, pass_good_negative = 1 - tau >= tau_threshold, metric = tau))
	}
	attr(f, "ObjCall") <- sprintf("messinaSurvObj.Tau(tau_threshold = %s)", as.character(tau_threshold))
	attr(f, "PlotThresh") <- c(tau_threshold, 0.5, 1 - tau_threshold)
	attr(f, "PlotLimits") <- c(0, 1)
	f
}


#' @export
messinaSurvObj.RelTau = function(reltau_threshold)
{
	stopifnot(reltau_threshold > 0.5 && reltau_threshold < 1)

	f <- function(x, y)
	{
		sort_perm = order(y[,1])
		time = y[sort_perm,1]
		event = y[sort_perm,2]
		x = x[sort_perm]
		counts = messinaSurvConcordance(as.logical(x), time, event)

		reltau = counts$concordant / (counts$concordant + counts$discordant)

		return(list(pass_good_positive = reltau >= reltau_threshold, pass_good_negative = 1 - reltau >= reltau_threshold, metric = reltau))
	}
	attr(f, "ObjCall") <- sprintf("messinaSurvObj.RelTau(reltau_threshold = %s)", as.character(reltau_threshold))
	attr(f, "PlotThresh") <- c(reltau_threshold, 0.5, 1 - reltau_threshold)
	attr(f, "PlotLimits") <- c(0, 1)
	f
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


messinaSurvTrain = function(x1, y, min_group_frac, objective)
{
	# For each cutoff possible in the training set, determine if the
	# objective function is satisfied.
	x_sorted = rle(sort(x1))$values
	if (length(x_sorted) == 1)	{ x_cutoffs = x_sorted } else { x_cutoffs = (x_sorted[-1] + x_sorted[-length(x_sorted)])/2 }
	x_cutoffs = c(-Inf, x_cutoffs, Inf)
	x_cutoffs = x_cutoffs[!duplicated(x_cutoffs)]
	cutoff_obj = sapply(x_cutoffs, function(cutoff) {
		result = objective((x1 > cutoff)*1, y)
		c(result$pass_good_positive, result$pass_good_negative, result$metric) })

	# Also find which cutoffs satisfy the min_group_frac requirement.
	cutoff_frac = sapply(x_cutoffs, function(cutoff) mean(x1 > cutoff))
	cutoff_frac = pmin(cutoff_frac, 1 - cutoff_frac)
	cutoff_frac_good = cutoff_frac >= min_group_frac

	# Combine the two: objective function and min_group_frac, for each
	# possible classifier direction.
	cutoff_obj_good_positive = cutoff_obj[1,] & cutoff_frac_good
	cutoff_obj_good_negative = cutoff_obj[2,] & cutoff_frac_good
	
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
	
	result = list(threshold = best_threshold, margin = best_margin, direction = best_direction, cutoffs = x_cutoffs, objective_value = cutoff_obj[3,], min_group_frac = min_group_frac, objective = objective)
	return(result)
}


messinaSurvSingleXDraw = function(x1, y, min_group_frac, objective, n_test)
{
	test_indices = sample.int(length(x1), n_test)
	train_indices = setdiff(1:length(x1), test_indices)
	x_test = x1[test_indices]
	x_train = x1[train_indices]
	y_test = y[test_indices,]
	y_train = y[train_indices,]

	# Train a MessinaSurv classifier on the training set
	fit = messinaSurvTrain(x_train, y_train, min_group_frac, objective)
	
	# Assess the classifier on the test set
	if (is.na(fit$threshold))
	{
		test_obj = list(pass_good_positive = NA, pass_good_negative = NA, metric = NA)
	}
	else
	{
		test_class = 1*(x_test*fit$direction > fit$threshold*fit$direction)
		test_obj = objective(test_class, y_test)
	}
	return(test_obj)
}


messinaSurvSingleX = function(x1, y, min_group_frac, objective, n_boot, n_test)
{
	replicate(n_boot, messinaSurvSingleXDraw(x1, y, min_group_frac, objective, n_test), simplify = TRUE)
}


#' @importFrom plyr aaply llply
messinaSurvTrainOnSubset = function(x, y, min_group_frac, objective, subset, parallel, silent)
{
	fits = llply(1:nrow(x), 
		function(xi) 
		{ 
			if (subset[xi])
			{
				fit = messinaSurvTrain(x[xi,], y, min_group_frac, objective)
				return(fit)
			}
			else
			{
				return(NULL)
			}
		}, .parallel = parallel, .progress = ifelse(parallel || silent, "none", "time"))
	
	return(fits)
}

