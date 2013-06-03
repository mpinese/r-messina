/* errors.cpp
 * Error codes and messages.
 *
 * Copyright 2013 Mark Pinese
 *
 * Licensed under the Eclipse Public License 1.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://opensource.org/licenses/eclipse-1.0
 *
 * Changelog:
 * 20080513	Moved from types.h
 * 			Wrote setErrorCode, getLastErrorCode.
 * 20080722	Changed license from AFL 3.0 to CPL 1.0.
 * 20130603 Changed licence from CPL 1.0 to EPL 1.0.
 */

#include "errors.h"



static const char *messages[] = {
		"No error.",
		"Attempted to execute an unimplemented feature.  Please report this error and the conditions that generated it.",
		"Unknown error code.  Please report this error and the conditions that generated it.",
		"Invalid command line parameters.  See usage information.",
		"File open error.  Ensure that the filename and permissions are correct.",
		"Input file read error.  Ensure that the file is available for reading.",
		"Invalid input file format.  Check the file matches the required formatting.",
		"Data structure is being reinitialised.  Please report this error and the conditions that generated it.",
		"Memory allocation error.  Increase available RAM or reduce experiment complexity.",
		"Missing data encountered.  Ensure that no observations are missing or NA.",
		"Invalid classifier parameters.  Check that constraints fall within the permitted ranges.",
		"Attempt to use an uninitialised structure.  Please report this error and the conditions that generated it.",
		"Cached operation attempted, but cache has not been initialised.  Please report this error and the conditions that generated it.",
		"Classifier test attempted without prior training run.  Please report this error and the conditions that generated it.",
		"Data structure is invalid.  Please report this error and the conditions that generated it.",
		"Testing assertion failed.  Please report this error and the conditions that generated it.",
		"Senseless classifier parameters.  Ensure that sensitivity and specificity constraints result in an AUC > 0.5."
};



STATUS setErrorCode(STATUS error_code)
{
	static STATUS old_code = OK;
	STATUS last_code = old_code;
	old_code = error_code;
	return(last_code);
}


STATUS getLastErrorCode()
{
	STATUS last_code;
	last_code = setErrorCode(OK);
	setErrorCode(last_code);
	return(last_code);
}


const char *getErrorMsg(STATUS error_code)
{
	if (error_code == OK)
		return messages[0];
	else if (error_code == ERR_NOT_IMPLEMENTED)
		return messages[1];
	else if (error_code <= ERR_MAX_PLACEHOLDER)
		return messages[2];
	else
		return messages[2 - int(error_code)];
}

const char *getLastErrorMsg()
{
	return(getErrorMsg(getLastErrorCode()));
}
