/* Data.cpp
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
 * 20080414	Fixed an error in the asserts in processCell that was crashing
 * 			debug builds.
 * 20080722	Changed license from AFL 3.0 to CPL 1.0.
 * 20121010	Gutted to only include functionality required for an R interface.
 */

// Keep VS 2005 quiet about security issues.
#ifdef _MSC_VER
#define _CRT_SECURE_NO_DEPRECATE
#endif

#include <stdlib.h>
#include <string.h>

#include <assert.h>

#include "Data.h"
//#include "errors.h"



Data::Data()
{
	reset();
}


void Data::reset()
{
	m_init = false;
	m_exprs = NULL;
	m_classes = NULL;
	m_ngenes = 0;
	m_nsamps = 0;
	m_nclass1 = 0;
}


STATUS Data::allocData()
{
	if (isInit())
		return ERR_ALREADY_INIT;

	if ((m_exprs = new uint16_t[m_ngenes * m_nsamps]) == NULL)
		return ERR_MALLOC;
	if ((m_classes = new bool[m_nsamps]) == NULL)
	{
		delete[] m_exprs;
		return ERR_MALLOC;
	}
	
	return OK;
}


STATUS Data::destroyData()
{
	if (m_exprs != NULL)
		delete[] m_exprs;
	if (m_classes != NULL)
		delete[] m_classes;

	reset();	
	return OK;
}

/*
STATUS Data::printSummary() const
{
	printf("\nSummary of Data object at 0x%p:\n", this);
	
	if (!isInit())
	{
		printf("  Not initialized\n");
		return OK;
	}

	const uint8_t n_genes_print = 6, n_samps_print = 6;
	uint8_t i, j;
	printf("  Initialized, %ld genes x %ld samples\n", m_ngenes, m_nsamps);
	printf("  Data:");
	for (i = 0; i < (m_ngenes < n_genes_print ? m_ngenes : n_genes_print); i++)
	{
		printf("\n    ");
		for (j = 0; j < (m_nsamps < n_samps_print ? m_nsamps : n_samps_print); j++)
			printf("% 6d ", m_exprs[makeIndex(i, j)]);
		if (m_nsamps > n_samps_print)
			printf("...");
	}
	if (m_ngenes > n_genes_print)
	{
		printf("\n    ");
		for (j = 0; j < (m_nsamps < n_samps_print ? m_nsamps : n_samps_print) + 1; j++)
			printf("   ... ");
	}

	return OK;
}
*/

// Exprs is assumed to be in sample-major order.
STATUS Data::readMemory(int32_t ngenes, int32_t nsamps, const uint16_t *exprs, const bool *classes)
{
	if (isInit())
		return ERR_ALREADY_INIT;

	if (ngenes == 0 || nsamps == 0)
		return ERR_BAD_PARAM;
	
	m_ngenes = ngenes;
	m_nsamps = nsamps;
	
	STATUS err;
	if ((err = allocData()) != OK)
		return err;
	
	int32_t gene, sample;
	for (sample = 0; sample < nsamps; sample++)
	{
		m_classes[sample] = classes[sample];
		for (gene = 0; gene < ngenes; gene++)
			m_exprs[makeIndex(gene, sample)] = exprs[sample * ngenes + gene];
	}

	m_init = true;
	return OK;
}
