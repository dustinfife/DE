#'Generates \code{iterations} randomly generated models then fits them to the
#'dataset
#'
#'This function randomly generates \code{iterations} models, where the user
#'specifies the value of \code{iterations}. Each of these models are fit to the
#'dataset provided and the fit is assessed. For those models that exceed a
#'threshold, the model is outputted as a pdf graphic and/or csv RAM matrix.
#'
#'Model Conditioned Data Elasticity (DE) is a procedure proposed by Fife,
#'Rodgers, and Mendoza (2013). Given a correlation matrix, the procedure
#'randomly generates path analysis models. Those models that fit better than
#'the one proposed are flagged (either outputted as a path diagram via
#'graphviz, or outputted as a RAM matrix in csv form) for later inspection.
#'
#'The DE procedure also allows user restrictions. For example, suppose one
#'variable is Age. Obviously, Age should probably not be endogenous, so the
#'user can specify Age as an endogenous variable. That is done by creating a
#'matrix where the columns correspond to "From", "To", and "Include." For
#'example, to specify that A must cause B, one would insert in the first row of
#'the matrix c("A", "B", "1"). To specify that nothing can cause a variable
#'(i.e., to make a variable exogenous), one would leave the "From" column as
#'"". For example, the Age example would have c("", "Age", "0").
#'
#'@param data a correlation matrix containing the data to be fit
#'@param paths the number of paths for each randomly generated model. Can
#'either be a scalar, or a vector of length two, where the first element is the
#'minimum number of paths and the second element is the maximum. See examples.
#'@param location file location where the results are to be stored (as either a pdf or a RAM matrix)
#'@param var.names the names of the variables
#'@param restrictions a matrix containing restrictions on the DE procedure. See
#'details section.
#'@param arrows specifies what proportion of the paths should be bi-directional
#'(correlational)
#'@param actual.model a RAM matrix (containing columns from, to, arrows,
#'values) of the proposed model. If specified, either the point estimate of
#'RMSEA or its upper limit is used to identify confounding models.
#'@param actual.fit the RMSEA of the actual model. Serves as a threshold for
#'when a random model is outputted.
#'@param iterations Number of models randomly generated.
#'@param N the sample size that generated the correlation matrix.
#'@param graphviz logical. When a model performs better than the specified
#'threshold, should a graphic of the model be outputted? Requires prior
#'installation of graphviz.
#'@param output.ram logical. When a model performs better than the specified
#'threshold, should a csv of the model be outputted?
#'@param openmx logical. Should \code{OpenMx} be used to fit the model? If
#'FALSE, the \code{sem} package is used instead. The default is FALSE since
#'OpenMx cannot be installed from CRAN.
#'@param check.identical.models logical. Should all models be unique? Increases
#'computational intensity, but prevents duplicates in the DE procedure.
#'@param limit The number of times the DE procedure should try to find a unique
#'model before quitting early. Defaults at 100.
#'@param recursive Should paths be limited to recursive paths?
#'@param upper.limit logical. Should the upper limit of a 90\% confidence
#'interval for RMSEA be used to find confounding models?
#'@param ...  Other arguments passed to the DE.graphviz function.
#'@param fix.var Should the variances of the exogenous variables be fixed?
#'@param fitIndex specifies which fit index to use for comparison. One of the following, either "RMSEA" or "RMSR."
#'@param corr.exo Should the exogenous variables be automatically correlated?
#'@return \item{model.fits}{A Matrix where the first column is the iteration number and the second column is the fit index of that model.}
#'@return \item{better}{An integer indicating the number of models that fit better}
#'@return \item{fitIndex}{A string (either "RMSEA" or "RMSR") indicating which fit index was used.}
#'@author Dustin Fife
#'@seealso See Also \code{\link{DE}}, \code{\link{DE.fit}}, \code{\link{print.DE.dist}}
#'@references Fife, D.A., Rodgers, J.L., & Mendoza, J. L. (2013). Model
#'conditioned data elasticity in path analysis: Assessing the "confoundability"
#'of the data. Manuscript submitted for publication (under second revision).
#'@examples
#'
#'	#### load the albanese dataset
#'data(albanese)
#'
#'	#### set restrictions (make Age exogenous)
#'restrictions = matrix(c("", "Age", 0), nrow=1, byrow=TRUE)
#'
#'	#### make a RAM matrix of the albanese model
#'alb = matrix(c("Age", "External", 1, .2,
#'				"Age", "Reflexive", 1, .2,
#'				"Age", "CPM", 1, .2,
#'				"CPM", "Reflexive", 1, .2,
#'				"CPM", "Mental", 1, .2,
#'				"Age", "Mental", 1,  .2), nrow=6, byrow=TRUE)
#'alb = data.frame(alb); names(alb) = c("From", "To", "Arrows", "Values")	
#'	
#'	#### do DE
#'FF = DE.dist(data=albanese, paths=c(6,7), N=366, restrictions=restrictions, 
#'		location="deleteme", var.names=names(albanese), 
#'		actual.model=alb, iterations=20, arrows=0, openmx=TRUE, 
#'		upper.limit=TRUE, iterations=50)
#' @rdname DE.dist
#' @export DE.dist
#' @import OpenMx sem MBESS igraph MASS

DE.dist <-
function(data=NULL, paths, location="", var.names=NULL, restrictions=NULL, arrows=.2, actual.model=NULL, actual.fit=NULL, iterations=1, N=1000, graphviz=FALSE, output.ram=TRUE, openmx=FALSE, check.identical.models=TRUE, limit=100, upper.limit=TRUE, recursive=TRUE, fix.var=FALSE, fitIndex=c("RMSEA", "RMSR"), corr.exo=FALSE, ...){
	var.names = sort(var.names)
	
	require(MBESS, quietly=TRUE)
	require(igraph, quietly=TRUE)
	require(OpenMx, quietly=TRUE)	
	
	fitIndex = match.arg(fitIndex, c("RMSEA", "RMSR"))
	
	if (fitIndex=="RMSR" & upper.limit){
		stop("The upper limit of RMSR cannot be computed. Please either specify RMSEA as the fit index, or set upper.limit to FALSE.")
	}
	
	if (fitIndex=="RMSR" & !openmx){
		stop("Computing RMSR is only implemented withing openmx currently.")
	}
		
	##### fit user specified model (if supplied)
	if (!is.null(actual.model)){
		actual.model = data.frame(actual.model, stringsAsFactors=FALSE)
		actual.model[,1] = as.character(actual.model[,1])
		actual.model[,2] = as.character(actual.model[,2])
		actual.model[,3] = as.numeric(actual.model[,3])
		actual.model[,4] = as.numeric(actual.model[,4])		
		actual.model = actual.model[order(actual.model$From, actual.model$To),]				
		user.model = try(DE.fit(actual.model, dataset=data, N, openmx, fix.variances = fix.var))
		#print(summary(user.model))
		if (openmx){
			if (fitIndex=="RMSEA"){
				actual.fit = summary(user.model)$RMSEA
			} else {
				actual.fit = RMSR(data, model=user.model)
			}
			ests = summary(user.model)$parameters
			rws.t = which(ests$row != ests$col)
			actual.model[,4] = summary(user.model)$parameters$Estimate[rws.t]
		} else {
			actual.fit = summary(user.model)$RMSEA[1]
			ests = user.model$coeff
			rws.t = grep("par", names(ests))
			actual.model[,4] = ests[rws.t]			
		}
		fitSubsForPlot = actual.fit
		####### replace actual value with upper limit		
		if (openmx){
			if (upper.limit){
				actual.fit = ci.rmsea(summary(user.model)$RMSEA, df = summary(user.model)$degreesOfFreedom, 
						N=summary(user.model)$numObs, conf.level=.90)$Upper.Conf.Limit
				fitSubsForPlot = summary(user.model)$RMSEA		
			}
		} else {
			if (upper.limit){
				actual.fit = summary(user.model)$RMSEA[3]				
				fitSubsForPlot = summary(user.model)$RMSEA[1]				
			}

		}				
				
		#### output to graphviz 
		if (graphviz){
			lab = paste(fitIndex, ":", round(fitSubsForPlot, digits=4))
			DE.graphviz(model=actual.model, file=paste(location,"SubstantiveModel", sep=""), label=lab,...)
		}
		
		
		
	}
	

	##### install OpenMx if they haven't already done so. 	
	openMXinstalled = "OpenMx" %in% rownames(installed.packages())
	if (!openMXinstalled){
		ANSWER = readline("The DE package requires the OpenMx package. Do you want to install OpenMx? (Type 1 for yes and 2 for no) ")
		if (substr(ANSWER,1,1) == "2")
			stop("Quitting DE fitting algorithm. If you still wish to use DE, rerun and choose to install OpenMx.")
		else 
			cat("\nInstalling OpenMx\n")
			source('http://openmx.psyc.virginia.edu/getOpenMx.R')
			cat("\nInstallation Complete")
	}
	

	cat("\nRunning DE. This may take awhile.\n")
 
	##### warnings
	if (is.null(data)){
		stop("You must specify a correlation or a covariance matrix to fit.")
	}
	if (is.null(var.names)){
		stop("You must specify the names of the variables.")
	}	
	if (is.null(actual.fit)){
		stop("You must specify an RMSEA value or supply a RAM matrix \n for the model of interest")
	}
	
	##### create directory if it doesn't exist
	dir.create(location, showWarnings=FALSE)

	##### preallocate fit matrix and all models
	model.fits = matrix(nrow=iterations, ncol=2)
	errors=0
	better = 0

	##### prepare paths output
	if (length(paths)==1){
		paths = c(paths, paths)
	}

	models = matrix(nrow=iterations, ncol=1)
	for (i in 1:iterations){

		a = 0
		##### check for identical models
		if (check.identical.models){

				#### on first iteration, don't worry about duplicating a model
			if (i==1){
				
				#### check if it's recursive
				if (recursive){
					recur = FALSE
					while (!recur){
						r.model = DE(variable.names=var.names, paths=round(runif(1, paths[1], paths[2])), restrictions=restrictions, prop.arrows=arrows, corr.exogenous=corr.exo)
						recur = is.recur(r.model, var.names)
					}
				}
				
				r.model = r.model[order(r.model$From, r.model$To),]
				models[i,1] = paste(r.model$From, r.model$To, r.model$Arrows, collapse="--")
				
			} else {
				
			#### keep looking for a new model until one is found
			#### or until you've tried limit times
				there.already = TRUE
				while (there.already & a<limit){
					
						#### if it must be recursive, keep searching until you find a recursive model
					if (recursive){
						recur = FALSE
						while (!recur){
							r.model = DE(var.names, round(runif(1, paths[1], paths[2])), restrictions=restrictions, arrows, corr.exogenous=corr.exo)
							recur = is.recur(r.model, var.names)
						}
					}
					
					#### now determine if this model has been used before
					r.model = r.model[order(r.model$From, r.model$To),]
					str = paste(r.model$From, r.model$To, r.model$Arrows, collapse="--")
					there.already = str %in% models[,1]
					a = a+1
				}
				models[i,1] = str
				
			}
			
			##### break out of loop early
			if (a==limit){
				cat(paste("\nI'm having trouble finding models that haven't been used.\nI'm quitting early.\nIteration Number:", i))
				break
			}
		}



		#### fit it
		mxMod = try(DE.fit(returned.model=r.model, dataset=data, N, openmx, fix.variances= fix.var))
		if (mode(mxMod)!="list" & mode(mxMod) != "S4"){
			errors = errors+1
		} else {

			##### record the fit of the model
			model.fits[i,1] = i
			if (openmx){ 
					#### compute fit index of interest
				if (fitIndex=="RMSEA"){
					newfit = summary(mxMod)$RMSEA
				} else {
					newfit = RMSR(data, model= mxMod)
				}
				
				model.fits[i,2] = newfit
			} else {
				model.fits[i,2] =  sqrt(summary(mxMod)$chisq-summary(mxMod)$df)/sqrt(summary(mxMod)$df*(N-1))
			}
			#### if fit surpasses actual.fit, then output 
			if (!is.na(model.fits[i,2])){
			if (model.fits[i,2]<=actual.fit){
				better=better+1
				if (openmx){
					ests = summary(mxMod)$parameters
					rws.t = which(ests$row != ests$col)
					r.model[,4] = summary(mxMod)$parameters$Estimate[rws.t]
				} else {
					ests = mxMod$coeff
					rws.t = grep("par", names(ests))
					r.model[,4] = ests[rws.t]
				}
						
				##### output model to graphviz
				if (graphviz){
					tit = paste(fitIndex, ": ", round(model.fits[i,2], digits=4), sep="")
					DE.graphviz(model=r.model, file=paste(location,"/randomModel", i, sep=""), label=tit, ...)
				} 
				
				##### output model to csv file
				if (output.ram){
					dir.create(paste(location, "/RamMatrices", sep=""), showWarnings=FALSE)
					write.csv(r.model, paste(location, "/RamMatrices/Model", i, ".csv", sep=""), row.names=FALSE)
				}	
			}}				

		}	

		if (i/100 == round(i/100)){
			cat(paste("\nCompleted ", i, " of ", iterations, " iterations.", sep=""))
		}
	}  #### end iterations loop

	cat(paste("\nNumber of models that fit better:", better, "\n", sep=""))
	
	#### output fits to a file
	m = data.frame(model.fits)
	names(m) = c("Model Number", "RMSEA")
	write.csv(m, paste(location, "/FitDistribution.csv", sep=""), row.names=FALSE)
	model.fits = data.frame(model.fits)
	names(model.fits) = c("Iteration", "RMSEA")
	
	# ###### output histogram
	# par(mar=c(2.5, 2.5, 3, 1), mgp=c(1.5, .5, 0), family="Times", font.main=1)	
	# xl = range(model.fits[,2], na.rm=TRUE)
	# hist(model.fits[,2], breaks=(nrow(model.fits)/10), main="Distribution of RMSEAs", xlab="RMSEA")
	# abline(v=actual.fit, lwd=3)

	final = list(fits=na.omit(model.fits), better = better, fit.index=fitIndex)
	attr(final, "class") = c("DE.dist")
	final

}


#' @title Print DE Summary
#' @name Print DE Summary
#' @rdname print
#' @method print DE.dist
#' @S3method print DE.dist
print.DE.dist = function(x, ...){
	better =(x$better)
	cat(paste("There were", better, "models that fit better\n\n"))
	cat("Quantiles of the fit distribution:\n")
	print(summary(x$fits[,2]))
}

#' @title Plot DE Distribution
#' @name Plot DE Distribution
#' @rdname plot
#' @method plot DE.dist
#' @S3method plot DE.dist
plot.DE.dist = function(x, y, ...){
	mn = paste("Distribution of", x$fit.index)
	hist(x$fits[,2], xlab=x$fit.index, main=mn,...)
}


