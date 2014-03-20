/* Classifier.cpp
 * A class to implement the threshold classifier.  Performs
 * both training and testing.
 *
 * Copyright 2014 Mark Pinese
 *
 * Licensed under the Eclipse Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://opensource.org/licenses/eclipse-1.0
 *
 * Changelog:
 * 20070809	Commenced writing.
 * 20070810	Completed writing; performed preliminary testing.
 * 20070920	Minor epiphany (no toilet) led to conception of regionTrain.
 * 20070927	Finally nutted out the findFeasibleRegion algorithm, and coded 
 * 			an implementation.
 * 20071002	Fixed a silly omission in train that caused the program to not 
 * 			use a ZeroR classifier in the case of an empty FR.
 * 20080414	Added licence header.
 * 20080414	Initialised Classifier members to reasonable defaults in the
 * 			constructor, to avoid unsightly unintialised values appearing
 * 			in the output in some cases.
 * 20080506	Moved init code from 20080514 to the train function, due to 
 * 			re-use of the class making the constructor code useless.
 * 20080506	Added train(bool sorted), trainOnCache and setupFullTrainCache
 * 			to simplify training on the full data set.
 * 20080722	Changed licence from AFL 3.0 to CPL 1.0.
 * 20130603 Changed licence from CPL 1.0 to EPL 1.0.
 *          Changed Classifier::decide to use the internal R PRNG. 
 *          Changed asserts to Rcpp stops.
 */

#include "Classifier.h"

#include <Rcpp.h>


STATUS Classifier::init(float targ_sens, float targ_spec, const Data *data)
{
	if (targ_sens < 0 || targ_sens > 1 || targ_spec < 0 || targ_spec > 1)
		return ERR_BAD_PARAM;
	
	m_targ_sens = targ_sens;
	m_targ_spec = targ_spec;
	m_data = data;
	m_perf.set = false;
	if ((m_cache = new uint16_t[data->m_nsamps]) == NULL)
		return ERR_MALLOC;
	if ((m_trainexprs = new uint16_t[data->m_nsamps]) == NULL)
	{
		delete m_cache;
		return ERR_MALLOC;
	}
	if ((m_trainclasses = new bool[data->m_nsamps]) == NULL)
	{
		delete m_cache;
		delete m_trainexprs;
		return ERR_MALLOC;
	}

	m_cached = false;
	m_type = TYPE_UNKNOWN;

	m_init = true;

	return OK;
}


STATUS Classifier::destroy()
{
	if (m_init)
	{
		delete m_cache;
		delete m_trainexprs;
		delete m_trainclasses;
	}
	m_init = false;
	return OK;
}


STATUS Classifier::cacheGene(int32_t gene)
{
	if (!m_init)
		return ERR_NOT_INIT;

	int32_t i;
	for (i = 0; i < m_data->m_nsamps; i++)
		m_cache[i] = m_data->m_exprs[m_data->makeIndex(gene, i)];
	
	m_cached = true;
	return OK;
}


void Classifier::setupFullTrainCache(bool sorted)
{
	int32_t i;
	for (i = 0; i < m_data->m_nsamps; i++)
	{
		m_trainexprs[i] = m_cache[i];
		m_trainclasses[i] = m_data->m_classes[i];
	}

	if (!sorted)
		sortTrainCache(m_data->m_nsamps);
}


void Classifier::setupTrainCache(const int32_t *samples, int32_t n_samples, bool sorted)
{
	int32_t i;
	for (i = 0; i < n_samples; i++)
	{
		m_trainexprs[i] = m_cache[samples[i]];
		m_trainclasses[i] = m_data->m_classes[samples[i]];
	}

	if (!sorted)
		sortTrainCache(n_samples);
}


void Classifier::sortTrainCache(int32_t n_samples)
{
	// Insertion sort.
	int32_t i, j;
	uint16_t temp_expr;
	bool temp_class;
	for (i = 1; i < n_samples; i++)
	{
		temp_expr = m_trainexprs[i];
		temp_class = m_trainclasses[i];
		for (j = i; j >= 1 && m_trainexprs[j - 1] > temp_expr; j--)
		{
			m_trainexprs[j] = m_trainexprs[j - 1];
			m_trainclasses[j] = m_trainclasses[j - 1];
		}
		m_trainexprs[j] = temp_expr;
		m_trainclasses[j] = temp_class;
	}
}


inline bool Classifier::doesPerfPass(int32_t tp, int32_t fp, int32_t tn, int32_t fn, bool k_is_pos)
{
	float sens, spec;
	int32_t k_tp, k_tn;
	if (k_is_pos)
	{
		k_tp = tp;
		k_tn = tn;
	}
	else
	{
		k_tp = fn;
		k_tn = fp;
	}
	sens = float(k_tp) / (tp + fn);
	spec = float(k_tn) / (tn + fp);
	return (sens >= m_targ_sens && spec >= m_targ_spec);
}


void Classifier::findFeasibleRegion(int32_t n_samples, bool k_is_pos, int32_t& f0, int32_t& f1)
{
	int32_t tp, fp, tn, fn, np, nn, i, lf0, lf1;
	bool perf_passes, ifr = false;
	
	// Calculate initial performance.  We start with i = 0, therefore 
	// if k_is_pos == true (k = 1), then for i = 0 all samples will be
	// considered positive.  Conversely, if k_is_pos == false (k = -1),
	// then for i = 0 all samples will be considered negative.
	np = 0;
	nn = 0;
	for (i = 0; i < n_samples; i++)
	{
		if (m_trainclasses[i])
			np++;
		else
			nn++;
	}
	
	if (k_is_pos)
	{
		tp = np;
		fp = nn;
		tn = 0;
		fn = 0;
	}
	else
	{
		tp = 0;
		fp = 0;
		tn = nn;
		fn = np;
	}
	
	// Set f0 to its default value of no feasible region
	lf0 = n_samples;
	
	for (i = 0; i <= n_samples; i++)
	{
		// Check performance for the current cutoff.
		// NB: Note that the k_is_pos argument to doesPerfPass has been
		// overridden.  This is as this function keeps a correct tally
		// of tp and fp, whereas train (for which doesPerfPass was
		// originally written) keeps a tally assuming k = 1.  Therefore,
		// when evaluating performance, train requires a correction for
		// the k = -1 case.  No such correction is needed here, and so
		// we always set k_is_pos to true, which prevents doesPerfPass
		// from performing any corrections.
		perf_passes = doesPerfPass(tp, fp, tn, fn, true); 
		if (perf_passes && !ifr)
		{
			// We've just entered the feasible region.  If i = 0, then
			// the feasible region is unbounded low.  We set the first
			// sample *out* of the feasible region to be the bound; this
			// is useful for the margin calculations, and produces the
			// reasonable value of -1.
			lf0 = i - 1;
			ifr = true;
		}
		else if (!perf_passes && ifr)
		{
			// We've just left the feasible region.  fr1 is calculated
			// outside of the loop, to account for the i = n + 1 possibility.
			break;
		}
		
		// Bail if we're at the final iteration -- there's no point in updating
		// performance as we've already analysed all samples.
		// Note that a continue is required to ensure that i = n + 1 if the loop
		// terminates naturally; see later code which relies upon this behaviour.
		if (i == n_samples)
			continue;
		
		// Update performance
		if (k_is_pos)
		{
			if (m_trainclasses[i])
			{
				tp--;
				fn++;
			}
			else
			{
				fp--;
				tn++;
			}
		}
		else
		{
			if (m_trainclasses[i])
			{
				fn--;
				tp++;
			}
			else
			{
				tn--;
				fp++;
			}
		}
	}
	
	// Calculate f1.  Placement of the calculation here allows for a 
	// natural handling of the unbounded high feasibility region case,
	// as at the natural termination of the loop i = n + 1, therefore
	// f1 will equal n, indicating the unbounded condition (the maximum
	// allowable value that corresponds to a sample is n - 1).
	lf1 = i - 1;
	
	// At this point either: we've broken from the loop as 
	// we've left the feasible region, the loop has terminated
	// because we've never left the feasible region, or the
	// loop has terminated because we've never entered the
	// feasible region.  The cases can be separated as follows:
	// Region entered?		Region left?		Condition
	// YES					YES					0 <= f0 < f1 < n_samples
	// YES					NO					0 <= f0 < f1 = n_samples
	// NO					NO					n_samples = f0 = f1 = n_samples
	//
	// Thus f0 and f1 encode all required data:
	// FR bounded low, 		FR bounded high		0 <= f0 < f1 < n_samples
	// FR unbounded low, 	FR bounded high		-1 = f0 < f1 < n_samples
	// FR bounded low, 		FR unbounded high	0 <= f0 < f1 = n_samples
	// FR unbounded low,	FR unbounded high	-1 = f0 < f1 = n_samples
	// No FR									n_samples = f0 = f1 = n_samples
	// If the low bound exists, it is defined by f0, and if the high bound
	// exists it is defined by f1.
	f0 = lf0;
	f1 = lf1;
}


bool Classifier::makeUnboundedClassifierFromFR(int32_t n_samples, bool k_is_pos, int32_t f0, int32_t f1)
{
	// Attempts to construct an unbounded classifier based upon
	// the feasible region bounds.  If such a classifier exists,
	// sets *this accordingly and returns true, else returns
	// false.
	if (f0 == -1)
	{
		// At least unbounded low.  May be unbounded high as well,
		// but we don't really care -- as an unbounded low classifier
		// will do, we'll just return one.
		// For k = 1, an unbounded low classifier always returns 
		// positive; for k = -1, it always returns negative.
		m_type = TYPE_ONE_CLASS;
		m_posk = k_is_pos;
		return true;
	}
	
	if (f0 != n_samples && f1 == n_samples)
	{
		// Unbounded high, bounded low.  An unbounded high classifier
		// will do in this case.
		// For k = 1, an unbounded high classifier always returns
		// negative; for k = -1, it always returns positive.
		m_type = TYPE_ONE_CLASS;
		m_posk = !k_is_pos;
		return true;
	}
	
	return false;
}


STATUS Classifier::trainOnCache(int32_t n_samples)
{
	if (!m_init)
		return ERR_NOT_INIT;
	if (!m_cached)
		return ERR_NOT_CACHED;
	
	// Set members to reasonable default values 
	m_type = TYPE_UNKNOWN;
	m_ptrue = -1;
	m_cutoff = 0;
	m_margin = 0;
	m_cutoff = 0;
	
	int32_t pos_f0, pos_f1, neg_f0, neg_f1;
	
	// Find the feasible region for k = 1
	findFeasibleRegion(n_samples, true, pos_f0, pos_f1);
	
	// Check for unbounded solutions for k = 1.
	if (makeUnboundedClassifierFromFR(n_samples, true, pos_f0, pos_f1))
		return OK;
	
	// Find the feasible region for k = -1
	findFeasibleRegion(n_samples, false, neg_f0, neg_f1);
	
	// Check for unbounded solutions for k = -1.
	if (makeUnboundedClassifierFromFR(n_samples, false, neg_f0, neg_f1))
		return OK;
	
	// Check for neither k yielding a FR
	// Added 20071002 -- silly omission
	if (neg_f0 == n_samples && neg_f1 == n_samples && pos_f0 == n_samples && pos_f1 == n_samples)
	{
		int32_t num_positive = 0, i;
		for (i = 0; i < n_samples; i++)
		{
			if (m_trainclasses[i])
				num_positive++;
		}		
		m_type = TYPE_ZERO_R;
		m_ptrue = float(num_positive) / n_samples;
		return OK;
	}

	// Calculate the margins
	uint16_t pos_margin, neg_margin;
	int32_t best_f0, best_f1;
	bool best_k;
	bool neg_valid_fr, pos_valid_fr;
	neg_valid_fr = !(neg_f0 == n_samples && neg_f1 == n_samples);
	pos_valid_fr = !(pos_f0 == n_samples && pos_f1 == n_samples);

	if (!(m_trainexprs[pos_f1] >= m_trainexprs[pos_f0]))
		Rcpp::stop("Internal messina assertion failed (m_trainexprs[pos_f1] >= m_trainexprs[pos_f0]).  Please report this to the package maintainer.");
	if (!(m_trainexprs[neg_f1] >= m_trainexprs[neg_f0]))
		Rcpp::stop("Internal messina assertion failed (m_trainexprs[neg_f1] >= m_trainexprs[neg_f0]).  Please report this to the package maintainer.");
	if (!(neg_valid_fr || pos_valid_fr))
		Rcpp::stop("Internal messina assertion failed (neg_valid_fr || pos_valid_fr).  Please report this to the package maintainer.");
	if (!(pos_f1 > pos_f0 || !pos_valid_fr))
		Rcpp::stop("Internal messina assertion failed (pos_f1 > pos_f0 || !pos_valid_fr).  Please report this to the package maintainer.");
	if (!(neg_f1 > neg_f0 || !neg_valid_fr))
		Rcpp::stop("Internal messina assertion failed (neg_f1 > neg_f0 || !neg_valid_fr).  Please report this to the package maintainer.");

	pos_margin = m_trainexprs[pos_f1] - m_trainexprs[pos_f0];
	neg_margin = m_trainexprs[neg_f1] - m_trainexprs[neg_f0];

	// Select the k that corresponds to the highest margin
	if ((pos_margin >= neg_margin && pos_valid_fr && neg_valid_fr) || (pos_valid_fr && !neg_valid_fr))
	{
		best_k = true;
		best_f0 = pos_f0;
		best_f1 = pos_f1;
	}
	else
	{
		best_k = false;
		best_f0 = neg_f0;
		best_f1 = neg_f1;
	}
	
	// Update the classifier parameters accordingly
	m_type = TYPE_THRESHOLD;
	m_cutoff = (m_trainexprs[best_f1] + m_trainexprs[best_f0]) / 2;
	m_margin = m_trainexprs[best_f1] - m_trainexprs[best_f0];
	m_posk = best_k;	

	return OK;
	
	/* The concept underlying this training technique follows.
	 * As the threshold is varied with a fixed k, sensitivity and specificity will 
	 * monotonically change.  Their direction of change will be reversed, and 
	 * dependent upon the data and the value of k.  However, we can at every point 
	 * decide upon an appropriate k to optimise performance, using the criterion
	 * in train (this_k_is_positive = tp > fp).  For a reversed k, the corrected
	 * sensitivity and specificity are 1 minus their original values.
	 * 
	 * Given the above, we can produce for every threshold position a plot of
	 * sensitivity and specificity vs threshold, given an optimal k.  This plot
	 * can be combined with our criteria to define an allowable region of 
	 * thresholds.  This allowable region effectively defines all locations which
	 * the threshold may occupy whilst still satisfying the performance constraints.
	 * Given the feasible region, it is reasonable to select its centre as the
	 * final threshold, thus maximising the margin to the start of the disallowed
	 * space.
	 * 
	 * A few niggly special cases exist.  One is in which the acceptable region is
	 * unbounded, either to one side, or on both.  In the former case an infinite
	 * margin (ie. one-class classifier) is appropriate.  In the latter, there's
	 * really very little we can do -- the user's supplied us with silly constraints.
	 * Another niggly case in that in which the feasible region doesn't exist.  This
	 * is easily handled by specifying a zero-rule classifier.
	 * 
	 * Now I just have to figure out how to actually implement it.
	 */
}


STATUS Classifier::train(const int32_t *samples, int32_t n_samples, bool sorted)
{
	// Prepare a local array of the training data, in sorted order.
	setupTrainCache(samples, n_samples, sorted);

	return trainOnCache(n_samples);
}


STATUS Classifier::train(bool sorted)
{
	// Prepare a local array of the training data, in sorted order.
	setupFullTrainCache(sorted);
	
	return trainOnCache(m_data->m_nsamps);
}



STATUS Classifier::test(const int32_t *samples, int32_t n_samples)
{
	if (!m_init)
		return ERR_NOT_INIT;
	if (!m_cached)
		return ERR_NOT_CACHED;
	if (m_type == TYPE_UNKNOWN)
		return ERR_NOT_TRAINED;

	int32_t i, sample_index;
	bool truth, prediction;

	int32_t fn = 0, fp = 0, tn = 0, tp = 0;

	for (i = 0; i < n_samples; i++)
	{
		sample_index = samples[i];
		prediction = decide(m_cache[sample_index]);
		truth = m_data->m_classes[sample_index];
		if (truth)
		{
			if (prediction)
				tp++;
			else
				fn++;
		}
		else
		{
			if (prediction)
				fp++;
			else
				tn++;
		}
	}

	m_perf.tpr = float(tp) / n_samples;
	m_perf.tnr = float(tn) / n_samples;
	m_perf.fnr = float(fn) / n_samples;
	m_perf.fpr = float(fp) / n_samples;
	m_perf.set = true;

	return OK;
}


inline bool Classifier::decide(uint16_t expr_level)
{
	switch (m_type)
	{
	case TYPE_UNKNOWN:
		return false;
	case TYPE_THRESHOLD:
		if (m_posk)
			return expr_level > m_cutoff;
		else
			return expr_level < m_cutoff;
	case TYPE_ONE_CLASS:
		return m_posk;
	case TYPE_ZERO_R:
		return Rcpp::runif(1, 0, 1)[0] < m_ptrue;
	}

	Rcpp::stop("Internal messina assertion failed: Classifier::decide fell through.  Please report this to the package maintainer.");
	return false;
}


void Classifier::updatePerformance(Perf& sum, Perf& sum_sq) const
{
	sum.fnr += m_perf.fnr;
	sum.tpr += m_perf.tpr;
	sum.fpr += m_perf.fpr;
	sum.tnr += m_perf.tnr;
	sum_sq.fnr += m_perf.fnr * m_perf.fnr;
	sum_sq.tpr += m_perf.tpr * m_perf.tpr;
	sum_sq.fpr += m_perf.fpr * m_perf.fpr;
	sum_sq.tnr += m_perf.tnr * m_perf.tnr;
}


STATUS Classifier::fillResults(Result &res)
{
	if (!m_init)
		return ERR_NOT_INIT;
	if (m_type == TYPE_UNKNOWN)
		return ERR_NOT_TRAINED;
	
	res.class_posk = m_posk;
	res.class_threshold = m_cutoff;
	res.class_ptrue = m_ptrue;
	res.class_type = m_type;
	res.class_margin = m_margin;
	
	return OK;
}


// For testing use only.
STATUS Classifier::testAssert(CLASSIFIER_TYPE type, bool k, uint16_t cutoff, uint16_t margin, float ptrue)
{
	if (m_type != type)
		return ERR_ASSERT_FAIL;
	
	switch (m_type)
	{
	case TYPE_THRESHOLD:
		if (cutoff != m_cutoff || k != m_posk || margin != m_margin)
			return ERR_ASSERT_FAIL;
		break;
	case TYPE_ZERO_R:
		if (ptrue != m_ptrue)
			return ERR_ASSERT_FAIL;
		break;
	case TYPE_ONE_CLASS:
		if (k != m_posk)
			return ERR_ASSERT_FAIL;
		break;
	case TYPE_UNKNOWN:
	default:
		break;
	}
	
	return OK;
}
