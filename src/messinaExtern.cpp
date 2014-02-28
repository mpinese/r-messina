/* messinaExtern.cpp
 * R interface for the Messina algorithm.
 *
 * Copyright 2014 Mark Pinese
 *
 * Licensed under the Eclipse Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://opensource.org/licenses/eclipse-1.0
 *
 * Changelog:
 * 20121010 Wrote.
 * 20130603 Placed under the EPL licence.
 * 20140228 Added progress and silent options.
 */

#include <Rcpp.h>

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

#include "types.h"
#include "Data.h"
#include "Classifier.h"
#include "crossval.h"
#include "errors.h"


using namespace Rcpp;
using namespace std;


// R interface declaration ///////////////////////////////////////////////////////////////////
extern "C" 
{
	SEXP messinaCextern(SEXP Rx, SEXP Rcls, SEXP Rn_boot, SEXP Rn_train, SEXP Rminsens, SEXP Rminspec, SEXP Rseed, SEXP Rprogress, SEXP Rsilent);
}

// Internal function declarations ////////////////////////////////////////////////////////////
STATUS convertRMatrix2Data(NumericMatrix &x, LogicalVector &cls, Data &data);
SEXP convertResults2R(Result *results, uint32_t n_results);


// Function definitions //////////////////////////////////////////////////////////////////////

/*	Launch point for Messina C code called from R.
	Rx			NumericMatrix	Data matrix, as integers in [0, 65535].  Probes / genes are in rows, samples in columns.
	Rcls		LogicalVector	Class membership vector, true/false for each sample.
	Rn_boot		Integer			Number of bootstrap iterations.
	Rn_train	Integer			Number of samples to draw for training in each bootstrap iteration.
	Rminsens	Float			Minimum classifier sensitivity.
	Rminspec	Float			Minimum classifier specificity.
	Rseed		Integer			PRNG seed.
	Rprogress	Logical			Display progress bar?
	Rsilent		Logical			Be completely silent (except for errors or warnings)?

	The Messina code populates an array of type Result (defined in Classifier.h), which then
	is converted back to SEXP objects for return to R.
*/
SEXP messinaCextern(SEXP Rx, SEXP Rcls, SEXP Rn_boot, SEXP Rn_train, SEXP Rminsens, SEXP Rminspec, SEXP Rseed, SEXP Rprogress, SEXP Rsilent)
{
	BEGIN_RCPP

	NumericMatrix x(Rx);
	LogicalVector cls(Rcls);
	
	uint32_t seed = as<uint32_t>(Rseed);
	uint32_t n_boot = as<uint32_t>(Rn_boot);
	uint32_t n_train = as<uint32_t>(Rn_train);
	float minsens = as<float>(Rminsens);
	float minspec = as<float>(Rminspec);	
	STATUS err;
	string errmsg;
	SEXP retval;
	bool progress = as<bool>(Rprogress);
	bool silent = as<bool>(Rsilent);
	
	RNGScope scope;	
	
	Data data;
	Result *results;
	Classifier classifier;
	FILE *out;
	
	if ((err = convertRMatrix2Data(x, cls, data)) != OK)
	{
		errmsg.assign(getErrorMsg(err));
		return wrap<string>(errmsg);
	}
	
	if ((err = classifier.init(minsens, minspec, &data)) != OK)
	{
		errmsg.assign(getErrorMsg(err));
		return wrap<string>(errmsg);
	}
	
	if ((results = new Result[data.getNGenes()]) == NULL)
	{
		errmsg.assign(getErrorMsg(ERR_MALLOC));
		return wrap<string>(errmsg);
	}
	
	if ((err = CrossVal::cv(n_train, n_boot, classifier, results, progress, silent)) != OK)
	{
		delete results;
		errmsg.assign(getErrorMsg(err));
		return wrap<string>(errmsg);
	}
	
	retval = convertResults2R(results, data.getNGenes());
	delete results;
	return retval;
	
	END_RCPP
}


/*	Convert an SEXP encoding an R matrix (passed as a reference to an Rcpp NumericMatrix)
	to the internal representation used by Messina, a Data object.  Coerces all entries of
	x to uint16_t type (bounds [0, 65535]), and higher R code must ensure that this bound is
	respected.
*/
STATUS convertRMatrix2Data(NumericMatrix &x, LogicalVector &cls, Data &data)
{
	R_len_t n_probes = x.nrow();
	R_len_t n_samples = x.ncol();
	R_len_t probe_i, sample_i;
	STATUS err;
	
	// Use the readMemory function of Data.  Must therefore have the entries of x
	// sitting in memory to be read.
	uint16_t *exprs_array;
	bool *class_array;
	
	// Allocate an array to hold the expression data
	if ((exprs_array = new uint16_t[n_probes * n_samples]) == NULL)
		return ERR_MALLOC;

	// Allocate an array to hold the class membership data
	if ((class_array = new bool[n_samples]) == NULL)
	{
		delete exprs_array;
		return ERR_MALLOC;
	}
	
	// Copy the data to their arrays
	for (sample_i = 0; sample_i < n_samples; sample_i++)
	{
		class_array[sample_i] = cls[sample_i];
		for (probe_i = 0; probe_i < n_probes; probe_i++)
		{
			exprs_array[probe_i + n_probes * sample_i] = static_cast<uint16_t>(x(probe_i, sample_i));	// Sample-major order
//			exprs_array[probe_i + n_probes * sample_i] = as<uint16_t>(x(probe_i, sample_i));	// Sample-major order
		}
	}
	
	// Initialise Data with the arrays
	err = data.readMemory(n_probes, n_samples, exprs_array, class_array);
	
	delete exprs_array;
	delete class_array;
	return err;
}


SEXP convertResults2R(Result *results, uint32_t n_results)
{
	IntegerMatrix d1(n_results, 3);		// class_type, class_threshold, class_margin
	NumericMatrix d2(n_results, 10);	// class_ptrue, p_successful, perf_mean(tpr,fpr,tnr,fnr), perf_var(tpr,fpr,tnr,fnr)
	LogicalVector d3(n_results);		// class_posk
	
	R_len_t i;
	Result *this_result;
	
	for (i = 0; i < n_results; i++)
	{
		this_result = results + i;
		
		d1(i, 0) = this_result->class_type;
		d1(i, 1) = this_result->class_threshold;
		d1(i, 2) = this_result->class_margin;
		
		d2(i, 0) = this_result->class_ptrue;
		d2(i, 1) = this_result->p_successful;
		d2(i, 2) = this_result->mean.tpr;
		d2(i, 3) = this_result->mean.fpr;
		d2(i, 4) = this_result->mean.tnr;
		d2(i, 5) = this_result->mean.fnr;
		d2(i, 6) = this_result->var.tpr;
		d2(i, 7) = this_result->var.fpr;
		d2(i, 8) = this_result->var.tnr;
		d2(i, 9) = this_result->var.fnr;
		
		d3[i] = this_result->class_posk;
	}
	
	return List::create(Named("d1") = d1, Named("d2") = d2, Named("d3") = d3);
}