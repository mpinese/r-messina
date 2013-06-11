messinaExtern <- function(Rx, Rcls, Rbootiter, Rn_train, Rminsens, Rminspec, Rseed)
{
	.Call("messinaCextern", Rx, Rcls, Rbootiter, Rn_train, Rminsens, Rminspec, Rseed, PACKAGE = "messina")
}


#setClass("MessinaResult", representation(a = "character", b = "numeric"))


messina = function(x, y, min_sens, min_spec, f_train = 0.9, n_boot = 50, seed = NULL)
{
	# x is an ExpressionSet or an object that can be converted to a matrix by as.matrix
	# y is a binary vector of class membership information
	
#	if (class(x) == "ExpressionSet")
#	{
#		features = featureNames(x)
#		x = exprs(x)
#	}
#	else
#	{
		x = as.matrix(x)
		features = rownames(x)
#	}

	# Calculate the number of training samples in each round.
	n_train = round(ncol(x) * f_train)
	n_train = min(max(2, n_train), ncol(x) - 1)
	
	# Basic error checking.
	target_AUC = (min_spec + min_sens) / 2
	if (!is.logical(y))					stop("Error: Class membership indicator y must be a vector of logical (TRUE / FALSE) values.")
	if (any(is.na(x)))					stop("Error: Missing values in x.  Messina cannot handle missing data at this time.")
	if (any(is.na(y)))					stop("Error: Missing values in y.  Messina cannot handle missing class membership information.")
	if (!is.numeric(min_sens))			stop(sprintf("Error: Invalid value for minimum acceptable sensitivity min_sens: %s", str(min_sens)))
	if (!is.numeric(min_spec))			stop(sprintf("Error: Invalid value for minimum acceptable specificity min_spec: %s", str(min_spec)))
	if (!is.numeric(f_train))			stop(sprintf("Error: Invalid value for training sample fraction f_train: %s", str(f_train)))
	if (!is.numeric(n_boot))			stop(sprintf("Error: Invalid value for number of bootstrap iterations n_boot: %s", str(n_boot)))
	if (ncol(x) < 3)					stop(sprintf("Error: Too few samples (%d).  Messina requires at least 3 samples.", ncol(x)))
	if (sum(y) == 0 || sum(!y) == 0)	stop(sprintf("Error: All samples are of the same class.  Messina requires samples in each class to operate."))
	if (ncol(x) != length(y))			stop(sprintf("Error: Inconsistent input data sizes.  Number of rows of data x (%d) does not equal length of class labels y (%d).", nrow(x), length(y)))
	if (min_sens <= 0 || min_sens > 1)	stop(sprintf("Error: Invalid value for minimum acceptable sensitivity parameter (%.3f).  Minimum acceptable sensitivity must be in (0, 1].", min_sens))
	if (min_spec <= 0 || min_spec > 1)	stop(sprintf("Error: Invalid value for minimum acceptable specificity parameter (%.3f).  Minimum acceptable specificity must be in (0, 1].", min_spec))
	if (f_train <= 0 || f_train >= 1)	stop(sprintf("Error: Invalid value for training sample fraction f_train (%.3f).  Training sample fraction must be in (0, 1), and preferably >= 0.5.", f_train))
	if (n_boot <= 0)					stop(sprintf("Error: Invalid value for number of bootstrap iterations (%d).  Number of bootstrap iterations must be >= 1 (and preferably >= 50).", f_train))
	if (target_AUC <= 0.5)				stop(sprintf("Error: Senseless parameters.  Supplied parameters (min_sens = %.3f, min_spec = %.3f) imply a target classifier AUC <= 0.5; all genes will satisfy these constraints.  Increase parameter values to make the classifiers more stringent.", min_sens, min_spec))

	# Sanity checks
	if (sum(y) < 5)						warning(sprintf("Warning: Low number of samples in positive class (%d).  Messina's performance may be unreliable for fewer than 5 samples per class.", sum(y)))
	if (sum(!y) < 5)					warning(sprintf("Warning: Low number of samples in negative class (%d).  Messina's performance may be unreliable for fewer than 5 samples per class.", sum(!y)))
	if (target_AUC < 0.7)				warning(sprintf("Warning: Very low stringency parameters.  Supplied parameters (min_sens = %.3f, min_spec = %.3f) imply a low target classifier AUC (AUC = %.3f); with these settings Messina's stringency will be particularly low.  Continuing anyway, but consider increasing parameter values to make the classifiers more stringent.", min_sens, min_spec, target_AUC))
	if (target_AUC >= 0.9)				warning(sprintf("Warning: Very high stringency parameters.  Supplied parameters (min_sens = %.3f, min_spec = %.3f) imply a high target classifier AUC (AUC = %.3f); with these settings Messina's stringency will be particularly high.  Continuing anyway, but consider using a conventional T-test or similar, which in this case may be more powerful than Messina.", min_sens, min_spec, target_AUC))
	if (f_train < 0.5)					warning(sprintf("Warning: Low training fraction f_train (%.3f).  Are you absolutely sure you want it this low?  Continuing anyway, but consider increasing f_train to >= 0.5 (recommended value is 0.9).", f_train))
	if (n_boot < 50)					warning(sprintf("Warning: Low number of bootstrap iterations n_boot (%d).  Are you absolutely sure you want it this low?  Continuing anyway, but consider increasing n_boot to >= 50 (recommended value is 50, more is better).", n_boot))
	
	# Transform x so that each gene (in rows) is bounded in [0, 65535].
	# Keep note of the transformation parameters so it can be reversed.
	xmin = apply(x, 1, min)
	xmax = apply(x, 1, max)
	xtrans = floor((x - xmin) / (xmax - xmin) * 65535)
	xtrans[xtrans == 65536] = 65535		# Just in case of rounding errors.
	stopifnot(all(xtrans <= 65535) & all(xtrans >= 0))
	
	# Get a random seed if one wasn't supplied
	if (is.null(seed))	seed = runif(1, 1, 2^32 - 1)
	
	result = messinaExtern(xtrans, y, n_boot, n_train, min_sens, min_spec, seed)
	
	if (typeof(result) == "character")
	{
		# messina encountered an error, which is returned here as a string.
		stop(sprintf("Error: Internal error encountered in messina C++ code.  Error description from messina: \"%s\".", result))
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
	classifier_threshold = result$d1[,2] * (xmax - xmin) / 65535 + xmin
	classifier_margin = result$d1[,3] * (xmax - xmin) / 65535
	classifier_threshold[classifier_type != "Threshold"] = NA
	classifier_margin[classifier_type != "Threshold"] = NA
	classifier_ptrue = result$d2[,1]
	classifier_ptrue[classifier_type != "Random"] = NA
	classifier_psuccessful = result$d2[,2]
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
	
	result2 = list(	problem = "classification",
					parameters = list(x = x, y = y, features = features, min_sens = min_sens, min_spec = min_spec, f_train = f_train, n_boot = n_boot, seed = seed),
					classifier = list(type = classifier_type, threshold = classifier_threshold, ptrue = classifier_ptrue, posk = classifier_posk), 
					margin = classifier_margin, 
					psuccessful = classifier_psuccessful, 
					passed = specifications_passed,
					perf = list(mean = classifier_perf_mean, var = classifier_perf_var))
	class(result2) = "MessinaResult"
	
	return(result2)
}



#plot.MessinaResult = function(result, i = NULL, type = "ggplot2")
plot.MessinaResult = function(x, ...)
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
	
	if (x$problem == "classification")
	{
		if (is.null(i))
		{
			temp.margin = x$margin
			temp.margin[!x$passed] = NA
			i = which.max(temp.margin)
		}
		this_x = x$parameters$x[i,]
		this_order = order(this_x)
		this_x = this_x[this_order]
		this_y = x$parameters$y[this_order]
		this_threshold = x$classifier$threshold[i]
		this_margin = x$margin[i]
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
				ggtitle(x$parameters$features[i]) + coord_flip()
		}
		else
		{		
			barplot(this_x, col = c("green", "red")[this_y + 1], ylim = c(0, ymax))
			abline(h = c(this_threshold - this_margin/2, this_threshold, this_threshold + this_margin/2), lwd = 2, lty = c("dotted", "solid", "dotted"))
		}
	}
	else if (x$problem == "survival")
	{
		stop("Error: Plotting for messina problem type 'survival' is not yet implemented.")
	}
	else
	{
		stop(sprintf("Error: Unrecognised problem type '%s' -- cannot plot this Messina result.", x$problem))
	}
}


summary.MessinaResult = function(object, ...)
{
	toptable(object, ...)
}


toptable = function(result, maxn = 10, minmar = 1)
{
	if (result$problem == "classification")
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
		return(invisible(summary))
	}
	else if (result$problem == "survival")
	{
		stop("Error: toptable for messina problem type 'survival' is not yet implemented.")
	}
	else
	{
		stop(sprintf("Error: Unrecognised problem type '%s' -- cannot print this Messina result.", result$problem))
	}
}
