survivalMaxMarginFitHelperC <- function (Rsurv_t, Rsurv_e, Rtimes, Rgroup) 
{
	.Call("survivalMaxMarginFitHelperCextern", Rsurv_t, Rsurv_e, Rtimes, Rgroup, PACKAGE = "messina")
}


survivalMaxMarginFitHelper = function(x, y)
{
	perm = order(x)
	y = y[perm,]
	x = x[perm]
	
	group = upper.tri(diag(length(x)), diag = FALSE)[-length(x),]*1
	fracs = rowMeans(group)
	cutoffs = x[-1]
	times = sort(unique(y[,1][y[,2] == TRUE]))
	J = length(times)
	K = nrow(group)

	ret = survivalMaxMarginFitHelperC(y[,1], y[,2], times, group)
	O = ret[["O"]]
	O0 = ret[["O0"]]
	N = ret[["N"]]
	N0 = ret[["N0"]]

	O1 = O - O0
	N1 = N - N0
	f1 = N1/N
	E1 = O*f1
	V = E1*(1-f1)*(N-O)/(N-1)
	V[N == 1] = 0
	Z = colSums(O1 - E1) / sqrt(colSums(V))
	p = pnorm(Z)
	p = pmin(p, 1 - p)*2
	list(cutoff = cutoffs, P = p)
}


messinaSurv = function(x, y, sig = 0.05, minfrac = 0.1, method = "fast")
{
	library(survival)
	if (method == "fast")
	{
		result = survivalMaxMarginFitHelper(x, y)
		cutoffs = result$cutoff
		pvals = result$P
		pvals[is.na(pvals)] = 1
	}
	else if (method == "orig")
	{
		xvals = sort(unique(x))
		cutoffs = xvals[-1]		# A if < cutoff, B if >= cutoff.
		pvals = sapply(cutoffs, function(cutoff) pchisq(survdiff(y ~ x < cutoff)$chisq, 1, lower.tail = FALSE))
	}
	else
	{
		stop(sprintf("Invalid method for survivalMaxMarginFit: %s", method))
	}
	
	fracs = sapply(cutoffs, function(cutoff) min(mean(x < cutoff), 1 - mean(x < cutoff)))
	pval_ok = pvals < sig & fracs >= minfrac
	
	# Now for the tricky bit: find the widest part of contiguous pval_ok == TRUE,
	# where width is determined by cutoffs.  For now take it slow and simple.
	last_ok = pval_ok[1]
	last_cutoff = cutoffs[1]
	widest_run_width = -1
	widest_run_cutoff = NA
	for (i in 2:length(pval_ok))
	{
		if (pval_ok[i] == TRUE & last_ok == FALSE)
		{
			# The start of a valid run.
			last_cutoff = cutoffs[i]
			last_ok = TRUE
		}
		else if (pval_ok[i] == FALSE & last_ok == TRUE)
		{
			# The end of a valid run.
			run_width = cutoffs[i] - last_cutoff
			if (run_width > widest_run_width)
			{
				widest_run_width = run_width
				widest_run_cutoff = last_cutoff + run_width/2
			}
			last_ok = FALSE
		}
	}
	# Note: isn't set up to handle edge cases at this time (ie open
	# intervals of pval_ok == TRUE).  Not much of a problem as these
	# aren't of any prognostic use anyway, especially thanks to the
	# minsize constraint.
	
	if (widest_run_width == -1)
	{
		return(c(NA, NA))
	}
	else
	{
		return(c(widest_run_cutoff, widest_run_width))
	}
}
