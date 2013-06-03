/* Code to implement the WELL1024 PRNG.
 * 
 * Modified from the original available at:
 * http://www.iro.umontreal.ca/~panneton/WELLRNG.html
 * accessed on 9 August 2007.
 * 
 * Mark Pinese
 * 20070809	Adapted from public code.
 * 20070928	Added randf, randF, randI, randi
 * 
 * Original copyright notice follows:
 */
/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */

#ifndef WELL1024A_H_
#define WELL1024A_H_

#include "types.h"

#define __WELL_CF_FULL	2.32830643653869628906e-10				// 1/(2^32)
#define __WELL_CF_PART	2.3283064359965952029459655278022e-10	// 1/(2^32 + 1)

namespace well
{

void seed(uint32_t seed);
void seeda(uint32_t *init);
uint32_t rand_raw();
inline double randf()		{ return (rand_raw() * __WELL_CF_PART);	}		// Samples from [0, 1)
inline double randF()		{ return (rand_raw() * __WELL_CF_FULL);	}		// Samples from [0, 1]
inline uint32_t randi(uint32_t upper_bound)	{ return (uint32_t(upper_bound * randf())); }
inline uint32_t randI(uint32_t upper_bound)	{ return (uint32_t(upper_bound * randF())); }

}

#endif /*WELL1024A_H_*/
