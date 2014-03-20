/* types.h
 * typedefs to handle varied treatment of integers by different 
 * platforms and compilers.
 *
 * Copyright 2014 Mark Pinese
 *
 * Licensed under the Eclipse Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://opensource.org/licenses/eclipse-1.0
 *
 * Changelog:
 * 20070709	Wrote.
 * 20080414	Added licence header.
 * 20080513	Moved error codes to errors.h
 * 20080722	Changed license from AFL 3.0 to CPL 1.0.
 * 20121010	Changed GNUC types to use stdint.h values.
 * 20130603 Changed licence from CPL 1.0 to EPL 1.0.
 */

#ifndef TYPES_H_
#define TYPES_H_

#if defined(__GNUC__)
	#include <stdint.h>
#elif defined(_MSC_VER)
	typedef __int64				int64_t;
	typedef __int32				int32_t;
	typedef __int16				int16_t;
	typedef __int8				int8_t;
	typedef unsigned __int64	uint64_t;
	typedef unsigned __int32	uint32_t;
	typedef unsigned __int16	uint16_t;
	typedef unsigned __int8		uint8_t;
#else
	#error "Unknown compiler or platform.  Please check and update types.h accordingly. "
#endif



#endif /*TYPES_H_*/
