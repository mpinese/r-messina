
#' @export
setMethod("show", signature = "MessinaResult", definition = function(object) 
{
	cat("An object of class ", class(object), "\n", sep = "")
	invisible(NULL)
})
