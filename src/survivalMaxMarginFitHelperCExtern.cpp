#include <Rcpp.h>

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;

extern "C" 
{
	SEXP survivalMaxMarginFitHelperCextern(SEXP Rsurv_t, SEXP Rsurv_e, SEXP Rtimes, SEXP Rgroup);
}


SEXP survivalMaxMarginFitHelperCextern(SEXP Rsurv_t, SEXP Rsurv_e, SEXP Rtimes, SEXP Rgroup)
{
	BEGIN_RCPP

	NumericVector surv_t(Rsurv_t);
	LogicalVector surv_e(Rsurv_e);
	NumericVector times(Rtimes);
	LogicalMatrix group(Rgroup);

	R_len_t n_surv = surv_t.size();
	R_len_t n_times = times.size();
	R_len_t n_groups = group.nrow();
	R_len_t time_i, group_i, surv_i;
	unsigned long accumO0, accumO, accumN0, accumN;

	NumericMatrix O0(n_times, n_groups);
	NumericMatrix O(n_times, n_groups);
	NumericMatrix N0(n_times, n_groups);
	NumericMatrix N(n_times, n_groups);
	for (time_i = 0; time_i < n_times; time_i++)
	{
		for (group_i = 0; group_i < n_groups; group_i++)
		{
			accumO0 = 0;
			accumO = 0;
			accumN0 = 0;
			accumN = 0;
			for (surv_i = 0; surv_i < n_surv; surv_i++)
			{
				accumO0 += group(group_i, surv_i) == false && surv_t(surv_i) == times(time_i) && surv_e(surv_i) == true;
				accumO  +=                                    surv_t(surv_i) == times(time_i) && surv_e(surv_i) == true;
				accumN0 += group(group_i, surv_i) == false && surv_t(surv_i) >= times(time_i);
				accumN  +=                                    surv_t(surv_i) >= times(time_i);
			}
			O0(time_i, group_i) = accumO0;
			O(time_i, group_i) = accumO;
			N0(time_i, group_i) = accumN0;
			N(time_i, group_i) = accumN;
		}
	}

	return (List::create(
				Named("O0") = O0,
				Named("O") = O,
				Named("N0") = N0,
				Named("N") = N));

	END_RCPP
}
