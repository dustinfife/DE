##' Compute Root Mean Square Residual (RMSR)
##'
##' @title Compute RMSR
##' @param data a correlation matrix
##' @param model an mxmodel object
##' @return the RMSR of the model
##' @author Dustin Fife
##' @export
RMSR = function(data, model){
	residuals <- cov2cor(data.matrix(data)) - cov2cor(model@objective@info$expCov)
	p = nrow(residuals)
	val = sqrt(sum(residuals^2)/(p*(p-1)))
	return(val)
}
