/* Code to implement the WELL1024 PRNG.
 * 
 * Modified from the original available at:
 * http://www.iro.umontreal.ca/~panneton/WELLRNG.html
 * accessed on 9 August 2007.
 * 
 * Mark Pinese
 * 20070809	Adapted from public code.
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

#include "well1024a.h"

#include "types.h"


namespace well
{

#define __WELL_MAT0POS(t,v) (v^(v>>t))
#define __WELL_MAT0NEG(t,v) (v^(v<<(-(t))))

#define __WELL_V0            __WELL_STATE[__WELL_state_i                   ]
#define __WELL_VM1           __WELL_STATE[(__WELL_state_i+3) & 0x0000001fU ]
#define __WELL_VM2           __WELL_STATE[(__WELL_state_i+24) & 0x0000001fU]
#define __WELL_VM3           __WELL_STATE[(__WELL_state_i+10) & 0x0000001fU]
#define __WELL_VRm1          __WELL_STATE[(__WELL_state_i+31) & 0x0000001fU]
#define __WELL_newV0         __WELL_STATE[(__WELL_state_i+31) & 0x0000001fU]
#define __WELL_newV1         __WELL_STATE[__WELL_state_i                   ]


static uint32_t __WELL_state_i = 0;
static uint32_t __WELL_STATE[32];
static uint32_t __WELL_z0, __WELL_z1, __WELL_z2;


void seed(uint32_t seed)
{
	// Just a bit of fun -- seed using the Collatz sequence
	// More serious users can employ seeda.
	int j;
	for (j = 0; j < 32; j++)
	{
		__WELL_STATE[j] = seed;
		if (seed % 2)
			seed = seed * 3 + 1;
		else
			seed = seed / 2;
	}
}


void seeda(uint32_t *init)
{
	int j;
	__WELL_state_i = 0;
    for (j = 0; j < 32; j++)
    	__WELL_STATE[j] = init[j];
}


uint32_t rand_raw()
{
	__WELL_z0 = __WELL_VRm1;
	__WELL_z1 = __WELL_V0 ^ __WELL_MAT0POS (8, __WELL_VM1);
	__WELL_z2 = __WELL_MAT0NEG (-19, __WELL_VM2) ^ __WELL_MAT0NEG(-14,__WELL_VM3);
	__WELL_newV1 = __WELL_z1 ^ __WELL_z2; 
	__WELL_newV0 = __WELL_MAT0NEG (-11,__WELL_z0) ^ __WELL_MAT0NEG(-7,__WELL_z1) ^ __WELL_MAT0NEG(-13,__WELL_z2);
	__WELL_state_i = (__WELL_state_i + 31) & 0x0000001fU;
	return (__WELL_STATE[__WELL_state_i]);
}


}
