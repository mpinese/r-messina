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
 * 201210
 */

#include <Rcpp.h>


using namespace Rcpp;
using namespace std;


// Internal function declarations ////////////////////////////////////////////////////////////
double coxScoreFunc(double theta, LogicalVector cls, NumericVector times, LogicalVector events);
double coxHessianFunc(double theta, LogicalVector cls, NumericVector times, LogicalVector events);


// Function definitions //////////////////////////////////////////////////////////////////////

/*	Fast implementation of a Cox model fit.  Arguments:
	cls			LogicalVector	Class membership vector, true/false for each sample.
	times		NumericVector	Observation time vector.
	events		LogicalVector	Event indicator vector (true = event, false = censored).
	tol 		double			Fitting terminates when abs(f'(theta)) <= tol
	maxiter		unsigned long  	or iter > maxiter.

	times should be sorted in increasing order.  Ties are currently not treated specially
	(ie "Breslow's method").

	Returns a single Numeric value, which is the fitted coefficient of the model
	coxph(Surv(times, events) ~ cls*1)
*/
// [[Rcpp::export]]
List messinaSurvCoxph(LogicalVector cls, NumericVector times, LogicalVector events, 
	double tol, unsigned long iter)
{
	double theta = 0;		// Starting coefficient estimate
	double deriv_at_theta;
	unsigned long iter;
	bool converged = false;

	for (iter = 0; iter < maxiter; iter++)
	{
		deriv_at_theta = coxScoreFunc(theta, cls, times, events);

		if (deriv_at_theta > -tol && deriv_at_theta < tol)
		{
			converged = true;
			break;
		}

		theta = theta - deriv_at_theta / coxHessianFunc(theta, cls, times, events);
	}

	return theta;
}


double coxScoreFunc(double theta, LogicalVector cls, NumericVector times, LogicalVector events)
{
	double score = 0;
	double inner_numer;
	unsigned long i, j;

	for (i = 0; i < cls.size(); i++)
	{
		if (events[i] == true)
		{
			if (cls[i] == true)
				score += 1;

			inner_numer = 0;
			for (j = i; j < cls.size(); j++)
			{
				if (cls[j] == true)
					inner_numer += theta;
			}
			score -= inner_numer / ()
		}
	}
}


double coxHessianFunc(double theta, LogicalVector cls, NumericVector times, LogicalVector events)
{

}
