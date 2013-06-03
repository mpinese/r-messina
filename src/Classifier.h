/* Classifier.h
 * A class to implement the threshold classifier.  Performs
 * both training and testing.
 *
 * Copyright 2008 Mark Pinese
 *
 * Licensed under the Common Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://www.opensource.org/licenses/cpl1.0.php
 *
 * Changelog:
 * 20070809	Commenced writing.
 * 20070910	Completed writing; performed preliminary testing.
 * 20070920	Minor epiphany (no toilet) led to writing of regionTrain.
 * 20080414	Added licence header.
 * 20080506	Added train(bool sorted), trainOnCache and setupFullTrainCache
 * 			to simplify training on the full data set.
 * 20080722	Changed license from AFL 3.0 to CPL 1.0.
 */

#ifndef _CLASSIFIER_H
#define _CLASSIFIER_H

#include "types.h"
#include "errors.h"
#include "Data.h"




enum CLASSIFIER_TYPE
{
	TYPE_UNKNOWN	= 0,
	TYPE_THRESHOLD	= 1,
	TYPE_ZERO_R		= 2,
	TYPE_ONE_CLASS	= 3
};


// A structure to hold performance metrics
struct Perf
{
	float tpr, fpr, tnr, fnr;
	bool set;
};


// A simple structure to store the main results of a classifier
// fit and its CV performance estimates.  The use of a smaller
// structure here avoids needing to allocate lots of bulky
// Classifier classes for all results.
struct Result
{
	CLASSIFIER_TYPE class_type;
	bool class_posk;
	float class_ptrue;
	uint16_t class_threshold, class_margin;
	
	float p_successful;
	Perf mean, var;
};


// Classifier is declared as a friend in Data, to substantially
// simplify its access to the underlying data.
class Classifier
{
	uint16_t m_cutoff;				// The threshold
	uint16_t m_margin;				// The associated margin
	bool m_posk;					// If m_posk == true, k = 1; else k = -1
	float m_ptrue;					// The probability of a true classification; used for the TYPE_ZERO_R classifier.
	CLASSIFIER_TYPE m_type;			// The type of this classifier.  Typically TYPE_THRESHOLD.
	float m_targ_sens, m_targ_spec;	// Target sensitivity and specificity values
	Perf m_perf;					// Performance figures

	const Data *m_data;				// The source data set
	uint16_t *m_cache;				// Gene expression value cache for rapid repeated access
	uint16_t *m_trainexprs;
	bool *m_trainclasses;

	bool m_init;					// Setup flag
	bool m_cached;

	bool decide(uint16_t expr_level);
	void setupFullTrainCache(bool sorted);
	void setupTrainCache(const int32_t *samples, int32_t n_samples, bool sorted);
	void sortTrainCache(int32_t n_samples);
	bool doesPerfPass(int32_t tp, int32_t fp, int32_t tn, int32_t fn, bool k_is_pos);
	void findFeasibleRegion(int32_t n_samples, bool k_is_pos, int32_t& f0, int32_t& f1);
	bool makeUnboundedClassifierFromFR(int32_t n_samples, bool k_is_pos, int32_t f0, int32_t f1);
	STATUS trainOnCache(int32_t n_samples);

public:
	Classifier(void)		{ m_init = false; };
	~Classifier(void)		{ destroy(); };

	STATUS init(float targ_sens, float targ_spec, const Data *data);
	STATUS destroy();
	STATUS cacheGene(int32_t gene);
	STATUS train(const int32_t *samples, int32_t n_samples, bool sorted);
	STATUS train(bool sorted);
	STATUS test(const int32_t *samples, int32_t n_samples);
	void updatePerformance(Perf& sum, Perf& sum_sq) const;
	
	const Data *getData() const			{ return m_data; };
	inline bool isInit() const			{ return m_init; };
	inline int32_t getNSamps() const		{ return (m_init ? m_data->m_nsamps : 0 ); };
	inline int32_t getNGenes() const		{ return (m_init ? m_data->m_ngenes : 0 ); };
	inline CLASSIFIER_TYPE getType() const	{ return m_type; };
	
	STATUS fillResults(Result &res);
	STATUS testAssert(CLASSIFIER_TYPE type, bool k, uint16_t cutoff, uint16_t margin, float ptrue);
};

#endif /*CLASSIFIER_H_*/
