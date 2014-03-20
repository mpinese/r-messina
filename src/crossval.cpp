/* crossval.cpp
 * Functions to implement cross-validation, in this case as
 * described by Nadeau & Bengio 2003.  Functions are wrapped
 * up in a namespace.
 *
 * Copyright 2014 Mark Pinese
 *
 * Licensed under the Eclipse Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://opensource.org/licenses/eclipse-1.0
 *
 * Changelog:
 * 20070812	Started writing.
 * 20070813	Completed writing.
 * 20080130	Added GUI reporting handling code.
 * 20080414	Added licence header.
 * 20080722	Changed license from AFL 3.0 to CPL 1.0.
 * 20130603 Changed licence from CPL 1.0 to EPL 1.0.
 *          Changed selectTestSet to use the internal R PRNG. 
 *          Changed asserts to Rcpp stops.
 * 20140131	Added progress bar code.
 * 20140228 Changed printf calls to Rprintf.
 *          Made progress bar optional with flag.
 */

#include <Rcpp.h>

#include "types.h"
#include "Classifier.h"
#include "Data.h"

#include "crossval.h"
#include "RInterrupt.h"

#include <stdio.h>


namespace CrossVal
{


inline void selectTestSet(bool *in_test, int32_t test_size, int32_t n_samples)
{
	int32_t test_done, new_point;

	if (!(test_size < n_samples))
		Rcpp::stop("Internal messina assertion failed (test_size < n_samples).  Please report this to the package maintainer.");
	
	for (new_point = 0; new_point < n_samples; in_test[new_point++] = false);
	
	for (test_done = 0; test_done < test_size; )
	{
		do
		{
			new_point = static_cast<int32_t>(Rcpp::floor(Rcpp::runif(1, 0, n_samples))[0]);
		}
		while (new_point == n_samples);		// According to the R docs this isn't required unless n_samples is 'small' compared to 0, but BSTS
		
		if (!in_test[new_point])
		{
			in_test[new_point] = true;
			test_done++;
		}
	}
}


STATUS gene_cv(int32_t train_size, uint16_t n_iters, Classifier& classifier, int32_t *train_indices, int32_t *test_indices, bool *in_test_set, Perf& mean, Perf& var, uint16_t& n_successful)
{
	uint16_t fold;
	int32_t sel_index, test_index = 0, train_index = 0, n_samps;
	STATUS err;
	Perf sum, sum_sq;
	float corr_factor;
	
	n_successful = 0;		// The total number of successful fits
	sum.tpr = 0;
	sum.tnr = 0;
	sum.fpr = 0;
	sum.fnr = 0;
	sum_sq.tpr = 0;
	sum_sq.tnr = 0;
	sum_sq.fpr = 0;
	sum_sq.fnr = 0;

	n_samps = classifier.getNSamps();

	// TODO: Consider making train_indices, test_indices, and in_test_set global within the namespace for a possible speed boost.
	for (fold = 0; fold < n_iters; fold++)
	{
		// Select new test points
		selectTestSet(in_test_set, n_samps - train_size, n_samps);

		// Convert the selection to indices
		test_index = 0;
		train_index = 0;
		for (sel_index = 0; sel_index < n_samps; sel_index++)
		{
			if (in_test_set[sel_index])
				test_indices[test_index++] = sel_index;
			else
				train_indices[train_index++] = sel_index;
		}

		// Train the classifier
		// TODO: Implement an early sort on the data, so that is_sorted can be passed as true, to avoid wasting lots of cycles in useless iterations of insertion sort.
		if ((err = classifier.train(train_indices, train_index, false)))
			return err;

		if (classifier.getType() == TYPE_THRESHOLD || classifier.getType() == TYPE_ONE_CLASS)
			n_successful++;
		
		// Test the classifier
		if ((err = classifier.test(test_indices, test_index)))
			return err;

		// Update the performance metrics
		// TODO: Will this call, which will be in the Classifier class, substantially slow things down?  Remember it's in the innermost CV loop.
		classifier.updatePerformance(sum, sum_sq);
	}

	// Use the performance sums to calculate mean and variance
	calcPerformanceStats(sum, sum_sq, mean, var, n_iters);

	// Calculate the Nadeau & Bengio correction factor
	// Here test_index is the test_size, and train_index is train_size
	corr_factor = 1 / float(n_iters) + (test_index / train_index);

	// Correct the variance estimates
	var.tpr *= corr_factor;
	var.tnr *= corr_factor;
	var.fpr *= corr_factor;
	var.fnr *= corr_factor;
		
	// Train on the full data set now
	if ((err = classifier.train(false)))
		return err;

	return OK;
}


inline void calcPerformanceStats(const Perf& sum, const Perf& sum_sq, Perf& mean, Perf& var, uint16_t n)
{
	mean.tpr = sum.tpr / n;
	mean.tnr = sum.tnr / n;
	mean.fpr = sum.fpr / n;
	mean.fnr = sum.fnr / n;

	var.tpr = (sum_sq.tpr*n - sum.tpr*sum.tpr) / (n*(n-1));
	var.tnr = (sum_sq.tnr*n - sum.tnr*sum.tnr) / (n*(n-1));
	var.fpr = (sum_sq.fpr*n - sum.fpr*sum.fpr) / (n*(n-1));
	var.fnr = (sum_sq.fnr*n - sum.fnr*sum.fnr) / (n*(n-1));
}


STATUS cv(int32_t train_size, uint16_t n_iters, Classifier& classifier, Result *results, bool progress, bool silent)
{
	// classifier.init has been called by the caller of this function.
	// results is assumed to be already allocated, and is a pointer to the start of an
	// array of length n_genes.
	
	bool *in_test_set;
	int32_t *train_indices, *test_indices;
	int32_t gene, n_genes, n_samps, test_size;
	STATUS err;
	Perf mean, var;
	uint16_t n_successful;
	const uint16_t progress_bar_width = 60;
	float progress_bar_frac_done;
	uint16_t progress_bar_n_done, progress_bar_i;

	if (!classifier.isInit())
		return ERR_NOT_INIT;
	
	n_genes = classifier.getNGenes();
	n_samps = classifier.getNSamps();
	
	if (train_size == 0 || train_size >= n_samps)
		return ERR_BAD_PARAM;
	test_size = n_samps - train_size;

	if ((train_indices = new int32_t[train_size]) == NULL)
		return ERR_MALLOC;

	if ((test_indices = new int32_t[test_size]) == NULL)
	{
		delete train_indices;
		return ERR_MALLOC;
	}

	if ((in_test_set = new bool[n_samps]) == NULL)
	{
		delete train_indices;
		delete test_indices;
		return ERR_MALLOC;
	}

	if (!silent)
		Rprintf("Performance bootstrapping...\n");

	for (gene = 0; gene < n_genes; gene++)
	{
		if (!silent && progress)
		{
			// Progress bar code
			if (gene % 100 == 0 || gene == n_genes - 1)
			{
				progress_bar_frac_done = float(gene + 1) / n_genes;
				progress_bar_n_done = uint16_t(progress_bar_frac_done * progress_bar_width);

				Rprintf("%3.0f%% [", progress_bar_frac_done * 100);
				for (progress_bar_i = 0; progress_bar_i < progress_bar_n_done; progress_bar_i++)
					Rprintf("=");

				for (; progress_bar_i < progress_bar_width; progress_bar_i++)
					Rprintf(" ");
				
				Rprintf("]\r");
			}
		}
		
		// User interrupt testing; run only every now and then.
		if (gene % 100 == 0)
		{
			if (RCheckInterrupt::checkInterrupt())
			{
				delete train_indices;
				delete test_indices;
				delete in_test_set;
				return ERR_ABORTED;
			}
		}
		
		if ((err = classifier.cacheGene(gene)) != OK)
		{
			delete train_indices;
			delete test_indices;
			delete in_test_set;
			return err;
		}
		
		if ((err = gene_cv(train_size, n_iters, classifier, train_indices, test_indices, in_test_set, mean, var, n_successful)) != OK)
		{
			delete train_indices;
			delete test_indices;
			delete in_test_set;
			return err;
		}
		
		if ((err = classifier.fillResults(results[gene])) != OK)
		{
			delete train_indices;
			delete test_indices;
			delete in_test_set;
			return err;
		}
		
		results[gene].mean = mean;
		results[gene].var = var;
		results[gene].p_successful = float(n_successful) / n_iters;
	}

	if (!silent)
		Rprintf("\n");

	delete train_indices;
	delete test_indices;
	delete in_test_set;

	return OK;
}

}
