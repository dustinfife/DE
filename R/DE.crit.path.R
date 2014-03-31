#'Find critical paths in a Path Analysis model.
#'
#'This function randomly generates models and determines which paths, when
#'removed, tend to reduce the fit of the model substantially.
#'
#'Fife, Rodgers, and Mendoza (2013) noted that DE distributions are frequently
#'bimodal and commented that a possible reason for this is that the poor
#'fitting distribution has one or more ``critical paths'' that are missing.
#'This function first generates \code{iterations} random models, then uses the
#'value for \code{cutpoint} to identify which paths exist in the right versus
#'left side of the cutpoint.
#'
#'@param data The correlation matrix to be fit.
#'@param variable.names A list of variables names if they're not specified in
#'the correlation matrix.
#'@param paths How many paths should be generated? Can be specified as a single
#'value or as a range. See examples.
#'@param prop.arrows What proportion of arrows should be correlational?
#'@param iterations Number of random models to be generated.
#'@param restrictions A matrix containing restrictions on the DE procedure. See
#'details section.
#'@param cutpoint An RMSEA value that differentiates between the modes of the
#'two distributions.
#'@return returns a matrix the specifies which paths are contained in the
#'"highMode" versus "lowMode" datasets
#'@export
#'@author Dustin Fife
#'@seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#'@references %% ~put references to the literature/web site here ~
#'@examples
#'
#'restrictions = matrix(c("", "Age", 0), nrow=1, byrow=TRUE)
#'data(albanese)
#'	
#'crit = DE.crit.path(data=albanese, paths=c(5,7), prop.arrows=0, iterations = 20, restrictions=restrictions, cutpoint=.4)
#'crit
#'
#' @import OpenMx sem MBESS igraph MASS
DE.crit.path <-
function(data, variable.names=names(data), paths, prop.arrows=.2, iterations=50, restrictions = NULL, cutpoint){
	
	##### create directories
	suppressWarnings(dir.create(paste(getwd(),"/critpath/", sep="")))
	suppressWarnings(dir.create(paste(getwd(),"/critpath/lowMode", sep="")))
	suppressWarnings(dir.create(paste(getwd(),"/critpath/highMode", sep="")))

	folders = getwd()
	
	###### iteratively fit models
	for (i in 1:iterations){
		model = DE(variable.names = names(cor), paths=paths, restrictions=restrictions, prop.arrows=prop.arrows)
		mxMod = try(DE.fit(model, dataset=cor, N=366, openmx=TRUE))
		if (mode(mxMod)=="list" | mode(mxMod) == "S4"){
			if (summary(mxMod)$RMSEA > cutpoint){
				write.csv(model, paste(folders, "/critpath/highMode", "/model", i, ".csv", sep=""), row.names=FALSE)
			} else {
				write.csv(model, paste(folders, "/critpath/lowMode", "/model", i, ".csv", sep=""), row.names=FALSE)
			}	
			
		}	
	}
	
	cat("Finished generating models. Now finding critical paths.\n")
	
	#### look at all possible connections
	possible = expand.grid(names(data), names(data), 1:2, stringsAsFactors=FALSE)
	possible[possible$Var1 == possible$Var2,] = NA
	possible = na.omit(possible)
	possible$numbHigh = 0
	possible$numbLow = 0
	allfilesHigh = list.files(paste(folders, "/critpath/highMode", sep=""))
	allfilesLow = list.files(paste(folders, "/critpath/lowMode", sep=""))
	
	###### warn if no variability
	if (length(allfilesHigh)<1 | length(allfilesLow)<1){
		stop("Based on your cutoff point, the distribution is not bimodal.\nTry adjusting your cutoff point.")
	}
	
	###### begin iteration
	for (i in 1:nrow(possible)){
		for (j in 1:length(allfilesHigh)){
			fl = read.csv(paste(folders, "/critpath/highMode/", allfilesHigh[j], sep=""), stringsAsFactors=FALSE)		
			pth = possible[i,]
			cond = which(pth$Var1 == fl$From & pth$Var2 == fl$To & pth$Var3 == fl$Arrows)
			if (length(cond)>0){
				possible[i,"numbHigh"] = possible[i,"numbHigh"] + 1
			}
		}
		for (k in 1:length(allfilesLow)){
			fl = read.csv(paste(folders, "/critpath/lowMode/", allfilesLow[k], sep=""), stringsAsFactors=FALSE)		
			pth = possible[i,]
			cond = which(pth$Var1 == fl$From & pth$Var2 == fl$To & pth$Var3 == fl$Arrows)
			if (length(cond)>0){
				possible[i,"numbLow"] = possible[i,"numbLow"] + 1
			}		
		}	
	}
	
	##### delete folders
	cmd2 = paste("rm -rf ", folders, "/critpath/highMode", sep="")
	try(system(cmd2))	
	cmd2 = paste("rm -rf ", folders, "/critpath/lowMode", sep="")
	try(system(cmd2))
	cmd2 = paste("rm -rf ", folders, "/critpath", sep="")
	try(system(cmd2))	


	##### find those rows where both are zero
	d = which(possible$numbHigh == 0 & possible$numbLow == 0)
	possible = possible[-d,]
	
	###### prep for output
	names(possible) = c("From", "To", "Arrows", "NumbHigh", "NumbLow")
	return(possible[order(possible$NumbHigh),])
}
