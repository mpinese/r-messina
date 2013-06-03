/* Data.h
 * A class to store transcriptional data, as well as class membership
 * for all samples in an experiment.  Also handles data IO.
 *
 * Copyright 2008 Mark Pinese
 *
 * Licensed under the Common Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://www.opensource.org/licenses/cpl1.0.php
 *
 * Changelog:
 * 20070809	Wrote and performed preliminary testing.
 * 20080414	Added licence header.
 * 20080722	Changed license from AFL 3.0 to CPL 1.0.
 * 20121010	Gutted to only include functionality required for an R interface.
 */

#ifndef DATA_H_
#define DATA_H_

#include <stdio.h>

#include "types.h"
#include "errors.h"


class Data
{
private:
	bool *m_classes;
	int32_t m_ngenes, m_nsamps, m_nclass1;
	bool m_init;
	uint16_t *m_exprs;

	void reset();
	int32_t makeIndex(int32_t gene, int32_t sample) const	{ return(gene * m_nsamps + sample); }; 
	STATUS allocData();
	STATUS destroyData();

	friend class Classifier;

public:
	Data();
	virtual ~Data()				{ destroyData(); };

	inline bool isInit() const	{ return m_init; };
	//STATUS printSummary() const;
	STATUS readMemory(int32_t ngenes, int32_t nsamps, const uint16_t *exprs, const bool *classes);

	inline int32_t getNGenes() const			{ return (m_init ? m_ngenes : 0); };
	inline int32_t getNSamples()	const		{ return (m_init ? m_nsamps : 0); };
	inline int32_t getNClass1Samples() const	{ return (m_init ? m_nclass1 : 0); };	
};

#endif /*DATA_H_*/
