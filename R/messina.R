#' Run the Messina algorithm to find features (eg. genes) that optimally distinguish
#' between two classes of samples, subject to minimum performance requirements.
#' 
#' TODO
#' 
#' @param x feature expression values, either supplied as an ExpressionSet, or as
#'   an object that can be converted to a matrix by as.matrix.  In the latter case,
#'   features should be in rows and samples in columns, with feature names taken
#'   from the rows of the object.
#' @param y a binary vector (TRUE/FALSE or 1/0) of class membership information for
#'   each sample in x.
#' @param min_sens the minimum acceptable sensitivity that a classifier separating
#'   the two groups of y must achieve.
#' @param min_spec the minimum acceptable specificity that a classifier separating
#'   the two groups of y must achieve.
#' @param f_train the fraction of samples to be used in the training splits of the
#'   bootstrap rounds.
#' @param n_boot the number of bootstrap rounds to use.
#' @param seed an optional random seed for the analysis.  If NULL, a random seed 
#'   derived from the current state of the PRNG is used.
#' @param progress display a progress bar tracking the computation?
#' @param silent be completely silent (except for error and warning messages)?
#' @return an object of class "MessinaClassResult" containing the results of the analysis.
#'
#' @export
#' @seealso \code{\link{MessinaClassResult-class}}
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @seealso \code{\link{messinaSurv}}
# @cite Pinese:2009
#' @author Mark Pinese \email{m.pinese@@garvan.org.au}
messina = function(x, y, min_sens, min_spec, f_train = 0.9, n_boot = 50, seed = NULL, progress = TRUE, silent = FALSE)
{
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

	# Calculate the number of training samples in each round.
	n_train = round(ncol(x) * f_train)
	n_train = min(max(2, n_train), ncol(x) - 1)
	
	progress = as.logical(progress)
	silent = as.logical(silent)
	if (is.na(progress))	{ progress = FALSE }
	if (is.na(silent))		{ silent = FALSE }
	
	# Basic error checking.
	target_AUC = (min_spec + min_sens) / 2
	if (!is.logical(y))					stop("Class membership indicator y must be a vector of logical (TRUE / FALSE) values.")
	if (any(is.na(x)))					stop("Missing values in x.  Messina cannot handle missing data at this time.")
	if (any(is.na(y)))					stop("Missing values in y.  Messina cannot handle missing class membership information.")
	if (!is.numeric(min_sens))			stop(sprintf("Invalid value for minimum acceptable sensitivity min_sens: %s", str(min_sens)))
	if (!is.numeric(min_spec))			stop(sprintf("Invalid value for minimum acceptable specificity min_spec: %s", str(min_spec)))
	if (!is.numeric(f_train))			stop(sprintf("Invalid value for training sample fraction f_train: %s", str(f_train)))
	if (!is.numeric(n_boot))			stop(sprintf("Invalid value for number of bootstrap iterations n_boot: %s", str(n_boot)))
	if (ncol(x) < 3)					stop(sprintf("Too few samples (%d).  Messina requires at least 3 samples.", ncol(x)))
	if (sum(y) == 0 || sum(!y) == 0)	stop(sprintf("All samples are of the same class.  Messina requires samples in each class to operate."))
	if (ncol(x) != length(y))			stop(sprintf("Inconsistent input data sizes.  Number of columns of data x (%d) does not equal length of class labels y (%d).", ncol(x), length(y)))
	if (min_sens <= 0 || min_sens > 1)	stop(sprintf("Invalid value for minimum acceptable sensitivity parameter (%.3f).  Minimum acceptable sensitivity must be in (0, 1].", min_sens))
	if (min_spec <= 0 || min_spec > 1)	stop(sprintf("Invalid value for minimum acceptable specificity parameter (%.3f).  Minimum acceptable specificity must be in (0, 1].", min_spec))
	if (f_train <= 0 || f_train >= 1)	stop(sprintf("Invalid value for training sample fraction f_train (%.3f).  Training sample fraction must be in (0, 1), and preferably >= 0.5.", f_train))
	if (n_boot <= 0)					stop(sprintf("Invalid value for number of bootstrap iterations (%d).  Number of bootstrap iterations must be >= 1 (and preferably >= 50).", f_train))
	if (target_AUC <= 0.5)				stop(sprintf("Senseless parameters.  Supplied parameters (min_sens = %.3f, min_spec = %.3f) imply a target classifier AUC <= 0.5; all genes will satisfy these constraints.  Increase parameter values to make the classifiers more stringent.", min_sens, min_spec))

	# Sanity checks
	if (sum(y) < 5)						warning(sprintf("Low number of samples in positive class (%d).  Messina's performance may be unreliable for fewer than 5 samples per class.", sum(y)))
	if (sum(!y) < 5)					warning(sprintf("Low number of samples in negative class (%d).  Messina's performance may be unreliable for fewer than 5 samples per class.", sum(!y)))
	if (target_AUC < 0.7)				warning(sprintf("Very low stringency parameters.  Supplied parameters (min_sens = %.3f, min_spec = %.3f) imply a low target classifier AUC (AUC = %.3f); with these settings Messina's stringency will be particularly low.  Continuing anyway, but consider increasing parameter values to make the classifiers more stringent.", min_sens, min_spec, target_AUC))
	if (target_AUC >= 0.9)				warning(sprintf("Very high stringency parameters.  Supplied parameters (min_sens = %.3f, min_spec = %.3f) imply a high target classifier AUC (AUC = %.3f); with these settings Messina's stringency will be particularly high.  Continuing anyway, but consider using a conventional T-test or similar, which in this case may be more powerful than Messina.", min_sens, min_spec, target_AUC))
	if (f_train < 0.5)					warning(sprintf("Low training fraction f_train (%.3f).  Are you absolutely sure you want it this low?  Continuing anyway, but consider increasing f_train to >= 0.5 (recommended value is 0.9).", f_train))
	if (n_boot < 50)					warning(sprintf("Low number of bootstrap iterations n_boot (%d).  Are you absolutely sure you want it this low?  Continuing anyway, but consider increasing n_boot to >= 50 (recommended value is 50, more is better).", n_boot))
	
	# Transform x so that each gene (in rows) is bounded in [0, 65535].
	# Keep note of the transformation parameters so it can be reversed.
	xmin = apply(x, 1, min)
	xmax = apply(x, 1, max)
	xdelta2 = xmax - xmin
	xdelta = xdelta2
	xdelta[xdelta == 0] = 1
	xtrans = floor((x - xmin) / xdelta * 65535)
	xtrans[xtrans == 65536] = 65535		# Just in case of rounding errors.
	stopifnot(all(xtrans <= 65535) & all(xtrans >= 0))
	
	# Call the external C functions for the actual calculation.
	result = messinaExtern(xtrans, y, n_boot, n_train, min_sens, min_spec, progress, silent)
	
	if (typeof(result) == "character")
	{
		# messina encountered an error, which is returned here as a string.
		stop(sprintf("Internal error encountered in messina C++ code.  Error description from C++ code: \"%s\".", result))
	}
	
	# Result is a list:
	#	Item		Dimensions		Data type	Description
	#   result$d1	Mx3 matrix		Integer		Columns are class_type, class_threshold, class_margin
	#	result$d2	Mx10 matrix		Numeric		Columns are class_ptrue, p_successful, mean tpr, mean fpr, mean tnr, mean fnr, var tpr, var fpr, var tnr, var fnr
	#	result$d3	M-vector		Logical		class_posk
	# M is the number of rows of x (genes / probes / whatever)
	# class_type takes the following values:
	#	0	UNKNOWN		Initialization value; shouldn't be encountered
	#	1	THRESHOLD	A conventional threshold classifier, one class if lower than / equal to the threshold, another if higher.
	#	2	ZERO_R		A marginal random classifier that returns class true if runif(1) < ptrue.
	#	3	ONE_CLASS	A deterministic single-class classifier that returns class posk.
	stopifnot(all(result$d1[,1] != 0))

	classifier_type = as.factor(c("Threshold", "Random", "OneClass")[result$d1[,1]])
	classifier_type[xdelta2 == 0] = "Random"
	classifier_threshold = result$d1[,2] * xdelta / 65535 + xmin
	classifier_margin = result$d1[,3] * xdelta / 65535
	classifier_threshold[classifier_type != "Threshold"] = NA
	classifier_margin[classifier_type != "Threshold"] = NA
	classifier_ptrue = result$d2[,1]
	classifier_ptrue[xdelta2 == 0] = mean(y)
	classifier_ptrue[classifier_type != "Random"] = NA
	classifier_psuccessful = result$d2[,2]
	classifier_psuccessful[xdelta2 == 0] = 0
	classifier_perf_mean = result$d2[,3:6]
	classifier_perf_var = result$d2[,7:10]
	classifier_posk = result$d3
	classifier_posk[classifier_type == "Random"] = NA
	
	colnames(classifier_perf_mean) = c("TPR", "FPR", "TNR", "FNR")
	colnames(classifier_perf_var) = c("TPR", "FPR", "TNR", "FNR")
	
	classifier_sens_mean = classifier_perf_mean[,"TPR"] / (classifier_perf_mean[,"TPR"] + classifier_perf_mean[,"FNR"])
	classifier_spec_mean = classifier_perf_mean[,"TNR"] / (classifier_perf_mean[,"TNR"] + classifier_perf_mean[,"FPR"])
	specifications_passed = classifier_sens_mean >= min_sens & classifier_spec_mean >= min_spec & classifier_type == "Threshold"
	classifier_perf_mean = cbind(classifier_perf_mean, Sensitivity = classifier_sens_mean, Specificity = classifier_spec_mean)
	
	params = .MessinaParameters(	x = x,
									y = y,
									features = features,
									samples = samples,
									perf_requirement = list(min_sensitivity = min_sens,
															min_specificity = min_spec),
									minimum_group_fraction = 0,
									training_fraction = f_train,
									num_bootstraps = n_boot,
									prng_seed = seed)
	
	perf_estimates = data.frame(mean_tpr = classifier_perf_mean[,"TPR"],
								mean_fpr = classifier_perf_mean[,"FPR"],
								mean_tnr = classifier_perf_mean[,"TNR"],
								mean_fnr = classifier_perf_mean[,"FNR"],
								var_tpr = classifier_perf_var[,"TPR"],
								var_fpr = classifier_perf_var[,"FPR"],
								var_tnr = classifier_perf_var[,"TNR"],
								var_fnr = classifier_perf_var[,"FNR"],
								mean_sens = classifier_sens_mean,
								mean_spec = classifier_spec_mean)
	
	fit_summary = data.frame(	passed = specifications_passed,
								type = classifier_type,
								threshold = classifier_threshold,
								posk = classifier_posk,
								margin = classifier_margin,
								ptrue = classifier_ptrue,
								psuccessful = classifier_psuccessful)
	
	fits = .MessinaFits(summary = fit_summary, objective_surfaces = list())
	
	result2 = .MessinaClassResult(	problem_type = "classification",
									parameters = params,
									perf_estimates = perf_estimates,
									fits = fits)
	
	
	return(result2)
}


messinaExtern <- function(Rx, Rcls, Rbootiter, Rn_train, Rminsens, Rminspec, Rprogress, Rsilent)
{
	.Call("messinaCextern", Rx, Rcls, Rbootiter, Rn_train, Rminsens, Rminspec, Rprogress, Rsilent, PACKAGE = "messina")
}
