##' Extract residuals of model
##'
##' Given a fitted model, \code{residuals} extracts the residual matrix for easier inspection
##' of the source of misfit.
##'	
##' @param object an \code{OpenMx} or \code{sem} model. 
##' @rdname residuals
##' @method residuals DE.fit
residuals.DE.fit = function(object,...){
	##### for sem package objects
	if (length(object$var.names)>0){
		residuals = residuals(object)
	} else {
		data = object$data@observed
		residuals <- cov2cor(data.matrix(data)) - cov2cor(object$objective@info$expCov)
	}
	return(residuals)
}