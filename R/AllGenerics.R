# AllGenerics.R:  Generic function definitions for Messina
# 
# Copyright 2014 Mark Pinese
#
# This file is distributed under the terms of the Eclipse Public 
# License v1.0, available at:
# https://www.eclipse.org/org/documents/epl-v10.html

#' @import graphics
if (!isGeneric("plot"))
{
	setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
}
