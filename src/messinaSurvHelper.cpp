/* messinaSurvHelper.cpp
 * Rcpp helper with fast survival functions.
 *
 * Copyright 2014-2015 Mark Pinese
 *
 * Licensed under the Eclipse Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://opensource.org/licenses/eclipse-1.0
 *
 * Changelog:
 * 20140730	Wrote first pass.
 * 20140805	Debugged coxphOptimizer.
 * 20150301	Removed coxphOptimizer and messinaSurvCoxph -- the fits
 *			were too slow, and occasionally unstable.  Replaced with
 *			messinaSurvLRT, which gives almost identical estimates, but
 *			is both much faster, and stable.
 */

#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <cmath>

using namespace Rcpp;
using namespace std;


// Internal function declarations ////////////////////////////////////////////////////////////


// Function definitions //////////////////////////////////////////////////////////////////////

/*	Fast implementation of a Log-rank test.  Arguments:
	cls			LogicalVector	Class membership vector, true/false for each sample.
	times		NumericVector	Observation time vector.
	events		LogicalVector	Event indicator vector (true = event, false = censored).

	times should be sorted in increasing order.

	Returns the test statistic, as a double.  The Z statistic is used; square this
	to get the standard 1 df ChiSq statistic.
*/
// [[Rcpp::export]]
double messinaSurvLRT(LogicalVector& cls, NumericVector& times, LogicalVector& events)
{
	unsigned long atrisk_1, atrisk, observed, observed_1;
	unsigned long n, i, j, k;
	double diff_sum, variance_sum, expected_1, variance;

	n = cls.size();

	atrisk = n;
	atrisk_1 = 0;
	for (i = 0; i < n; i++)
		atrisk_1 += cls[i] == false;

	diff_sum = 0;
	variance_sum = 0;
	for (i = 0; i < n; )
	{
		observed = 0;
		observed_1 = 0;
		for (j = i; j < n && times[j] == times[i]; j++)
		{
			observed += events[j] == true;
			observed_1 += events[j] == true && cls[j] == false;
		}

		expected_1 = double(observed) * double(atrisk_1) / double(atrisk);
		if (atrisk > 1)
			variance = expected_1 * (1 - double(atrisk_1) / double(atrisk)) * (double(atrisk) - double(observed)) / (double(atrisk) - 1);
		else
			variance = 0;

		diff_sum += double(observed_1) - expected_1;
		variance_sum += variance;

		atrisk -= (j - i);
		for (k = i; k < j; k++)
			atrisk_1 -= cls[k] == false;

		i = j;
	}

	return diff_sum / sqrt(variance_sum);
}


/*	Fast concordance calculation.  Arguments:
	cls			LogicalVector	Class membership vector, true/false for each sample.
	times		NumericVector	Observation time vector.
	events		LogicalVector	Event indicator vector (true = event, false = censored).

	times should be sorted in increasing order.

	Returns the a list, with members:
	"ties": count of comparisons with either tied times or class
	"concordant": count of concordant comparisons
	"discordant": count of discordant comparisons
*/
// [[Rcpp::export]]
List messinaSurvConcordance(LogicalVector& cls, NumericVector& times, LogicalVector& events)
{
	unsigned long ties = 0, concordant = 0, discordant = 0;
	unsigned long n, i, j;
	List ret;

	n = cls.size();

	for (i = 0; i < n; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			// Check for incomparability; skip if so.
			if ((times[i] < times[j] && events[i] == false) || (events[i] == false && events[j] == false))
				continue;

			if (cls[i] == cls[j] || times[i] == times[j])
				ties++;
			else if (cls[i] == true)
				concordant++;
			else
				discordant++;
		}
	}

	ret["ties"] = ties;
	ret["concordant"] = concordant;
	ret["discordant"] = discordant;

	return ret;
}
