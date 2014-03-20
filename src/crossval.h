/* crossval.h
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
 * 20140228 Added progress and silent options to cv.
 */

#ifndef _CROSSVAL_H
#define _CROSSVAL_H

#include "types.h"
#include "errors.h"
#include "Data.h"


namespace CrossVal
{


inline void calcPerformanceStats(const Perf& sum, const Perf& sum_sq, Perf& mean, Perf& var, uint16_t n);
inline void selectTestSet(bool *in_test, int32_t test_size, int32_t n_samples);
STATUS gene_cv(int32_t train_size, uint16_t n_iters, Classifier& classifier, int32_t *train_indices, int32_t *test_indices, bool *in_test_set, Perf& mean, Perf& var, uint16_t& n_successful);
STATUS cv(int32_t train_size, uint16_t n_iters, Classifier& classifier, Result *results, bool progress, bool silent);


}

#endif /*CROSSVAL_H*/
