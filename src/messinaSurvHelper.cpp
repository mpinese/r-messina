/* messinaSurvHelper.cpp
 * Rcpp helper with fast survival functions.
 *
 * Copyright 2014 Mark Pinese
 *
 * Licensed under the Eclipse Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://opensource.org/licenses/eclipse-1.0
 *
 * Changelog:
 * 20140730	Wrote first pass.
 * 20140805	Debugged coxphOptimizer.
 */

#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <cmath>

using namespace Rcpp;
using namespace std;


const double COX_DERIV2_MIN = 1.0e-6;		// The smallest absolute value of the Cox fit 2nd derivative permitted
const double COX_BETA_DELTA_MAX = 1.0e6;	// The largest absolute value of the Cox f'/f'' permitted.


// Internal function declarations ////////////////////////////////////////////////////////////
int coxphOptimizer(LogicalVector& cls, NumericVector& times, LogicalVector& events, double& beta, double tol, unsigned long maxiter, double sor_factor);


// Function definitions //////////////////////////////////////////////////////////////////////

/*	Fast implementation of a Cox model fit with a single binary independent variable.  Arguments:
	cls			LogicalVector	Class membership vector, true/false for each sample.
	times		NumericVector	Observation time vector.
	events		LogicalVector	Event indicator vector (true = event, false = censored).
	tol 		double			Fitting terminates when abs(f'(beta)) <= tol
	maxiter		unsigned long  	or iter > maxiter.

	times should be sorted in increasing order.  Ties are currently not treated specially
	(ie "Breslow's method").

	Returns a List, with members:
	  "beta"		the fitted coefficient of the model	coxph(Surv(times, events) ~ cls*1)
	  "converged"	an integer convergence flag.  0 => successful convergence.
*/
// [[Rcpp::export]]
List messinaSurvCoxph(LogicalVector& cls, NumericVector& times, LogicalVector& events, double tol, unsigned long maxiter)
{
	double beta = 0;	// Starting coefficient estimate
	int result;
	List ret;

	ret["converged"] = 0;

	if ((result = coxphOptimizer(cls, times, events, beta, tol, maxiter, 1.0)) == 0)
	{
		ret["beta"] = beta;
		return ret;
	}

	beta = 0;
	if ((result = coxphOptimizer(cls, times, events, beta, tol, maxiter, 0.5)) == 0)
	{
		ret["beta"] = beta;
		return ret;
	}

	RNGScope scope;
	beta = R::rnorm(0, 1);
	if ((result = coxphOptimizer(cls, times, events, beta, tol, maxiter, 0.5)) == 0)
	{
		ret["beta"] = beta;
		return ret;
	}

	ret["beta"] = NA_REAL;
	ret["converged"] = result;
	return ret;
}


/*	Fast implementation of a Log-rank test.  Arguments:
	cls			LogicalVector	Class membership vector, true/false for each sample.
	times		NumericVector	Observation time vector.
	events		LogicalVector	Event indicator vector (true = event, false = censored).

	times should be sorted in increasing order.

	Returns the test statistic, as a double.
*/
// [[Rcpp::export]]
double messinaSurvLRT(LogicalVector& cls, NumericVector& times, LogicalVector& events)
{
	unsigned long n1, n, i, o1, o, j, k;
	double diff_sum, var_sum, e1, var;

	n = cls.size();

	n1 = 0;
	for (i = 0; i < n; i++)
		n1 += cls[i] == false;

	diff_sum = 0;
	var_sum = 0;
	for (i = 0; i < n; )
	{
		for (j = i + 1; j < n && times[j] == times[i]; j++);

		o = j - i;
		o1 = 0;
		for (k = i; k < j; k++)
			o1 += cls[k] == false;

		e1 = double(o) / double(n) * n1;
		var = e1 * (1 - double(n1) / double(n))*(double(n) - double(o)) / (double(n) - 1);

		diff_sum += double(o1) - e1;
		var_sum += var;

		i = j;
	}

	return diff_sum / sqrt(var_sum);
}


int coxphOptimizer(LogicalVector& cls, NumericVector& times, LogicalVector& events, double& beta, double tol, unsigned long maxiter, double sor_factor)
{
	double deriv_at_beta, deriv2_at_beta, deriv_frac;
	unsigned long n_pos_events, n_pos, n;
	unsigned long i, iter, temp;

	// Get the number of observations in total
	n = cls.size();

	// Count the total number of events in group 1,
	// and the total size of group 1.
	n_pos_events = 0;
	n_pos = 0;
	for (i = 0; i < n; i++)
	{
		if (cls[i] == true)
		{
			n_pos++;
			if (events[i] == true)
				n_pos_events++;
		}
	}

	// Create a running count of the number in group 1
	// at and after time index i.
	vector<unsigned long> m_pos_i(n, n_pos);
	temp = 0;
	for (i = 0; i < n; i++)
	{
		m_pos_i[i] -= temp;
		if (cls[i] == true)
			temp++;
	}

	// Perform a Newton-Raphson search for the maximum likelihood.
	vector<double> deriv_temp(n);
	double exp_beta;
	double deriv_temp2;
	for (iter = 0; iter < maxiter; iter++)
	{
		// Calculate l'(beta)
		exp_beta = exp(beta);
		deriv_at_beta = double(n_pos_events);
		for (i = 0; i < n; i++)
		{
			if (events[i] == true)
			{
				deriv_temp2 = exp_beta * m_pos_i[i];
				deriv_temp[i] = deriv_temp2 / (deriv_temp2 + (n - i - m_pos_i[i]));
				deriv_at_beta -= deriv_temp[i];
			}
		}

		// Check for convergence
		if (deriv_at_beta > -tol && deriv_at_beta < tol)
			return 0;

		// Calculate l''(beta)
		deriv2_at_beta = 0;
		for (i = 0; i < n; i++)
		{
			if (events[i] == true)
				deriv2_at_beta -= deriv_temp[i]*(1 - deriv_temp[i]);
		}

		// Check for poor conditioning
		if (deriv2_at_beta > -COX_DERIV2_MIN && deriv2_at_beta < COX_DERIV2_MIN || !isnormal(deriv2_at_beta))
			return 2;
		deriv_frac = deriv_at_beta / deriv2_at_beta;
		if (deriv_frac < -COX_BETA_DELTA_MAX || deriv_frac > COX_BETA_DELTA_MAX || !isnormal(deriv_frac))
			return 3;

		// Update the search
		beta = (1-sor_factor)*beta + sor_factor*(beta - deriv_frac);
	}

	// Check for convergence after exactly maxiter iterations
	if (deriv_at_beta > -tol && deriv_at_beta < tol)
		return 0;

	return 1;
}
