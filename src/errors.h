/* errors.h
 * Error codes and messages.
 *
 * Copyright 2014 Mark Pinese
 *
 * Licensed under the Eclipse Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://opensource.org/licenses/eclipse-1.0
 *
 * Changelog:
 * 20080513	Moved from types.h
 * 20080722	Changed license from AFL 3.0 to CPL 1.0.
 * 20130603 Changed licence from CPL 1.0 to EPL 1.0.
 * 20130604 Added ERR_ABORTED value.
 */

#ifndef ERRORS_H_
#define ERRORS_H_


#define FP_ERROR_RETURN		-1000000000.0f


enum STATUS
{
	OK					= 0,
	ERR_BAD_CMDLINE		= -1,
	ERR_FILE_OPEN		= -2,
	ERR_FILE_READ		= -3,
	ERR_FILE_FORMAT		= -4,
	ERR_ALREADY_INIT	= -5,
	ERR_MALLOC			= -6,
	ERR_MISSING_DATA	= -7,
	ERR_BAD_PARAM		= -8,
	ERR_NOT_INIT		= -9,
	ERR_NOT_CACHED		= -10,
	ERR_NOT_TRAINED		= -11,
	ERR_BAD_DATASET		= -12,
	ERR_ASSERT_FAIL		= -13,
	ERR_BAD_CONSTRAINTS = -14,
	ERR_ABORTED			= -15,
	ERR_MAX_PLACEHOLDER = -16,
	ERR_NOT_IMPLEMENTED	= -100
};


const char *getErrorMsg(STATUS error_code);
STATUS setErrorCode(STATUS error_code);
STATUS getLastErrorCode();
const char *getLastErrorMsg();


#endif /*ERRORS_H_*/
