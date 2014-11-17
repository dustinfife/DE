#'Fit a RAM matrix to the data using path analysis.
#'
#'Given a single RAM matrix, this function fits the model and returns the
#'output from either \code{OpenMx} or \code{sem}.
#'
#'%% ~~ If necessary, more details than the description above ~~
#'
#'@param returned.model The model of interest, in RAM format. (A column of
#'"From" variables, "To" variables, the number of arrows, and the starting
#'value.)
#'@param dataset A correlation matrix.
#'@param N The sample sized used to estimate the correlation matrix.
#'@param openmx logical. Should \code{OpenMx} be used to fit the model? If
#'FALSE, the \code{sem} package is used instead. The default is FALSE since the
#'\code{sem} package is faster.
#'@param fix.variances logical. Should variances be fixed to one?
#'@param calc.se Should the standard errors be calculated? Defaults to no to speed up the DE algorithm. 
#'@return Depending on which program is used to fit the data, either an OpenMx
#'or an sem output.
#'@note %% ~~further notes~~
#'@author Dustin Fife
#'@seealso \code{\link{DE}}, \code{\link{DE.dist}}
#'@references %% ~put references to the literature/web site here ~
#'@keywords ~kwd1 ~kwd2
#'@examples
#'
#'data(albanese)
#'
#'	#### make restrictions
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
#'model = DE.fit(returned.model=alb, dataset=albanese, N = 366, openmx=TRUE, fix.variances=TRUE)
#'summary(model)
#'@export
#' @import OpenMx sem MBESS igraph MASS
DE.fit <-
function(returned.model, dataset, N=1000, openmx=FALSE, fix.variances=FALSE, calc.se=FALSE){
	b = returned.model

	##### convert to numeric and factor where needed
	b = data.frame(b)
	names(b) = c("From", "To", "Arrows", "Values")
	b$From = as.character(b$From); b$To = as.character(b$To); b$Arrows = as.numeric(as.character(b$Arrows)); b$Values = as.numeric(as.character(b$Values))
	
	##### strip illegal characters
	b$From = gsub("[^[:alnum:]]", "", b$From)
	b$To = gsub("[^[:alnum:]]", "", b$To)	
	names(dataset) = gsub("[^[:alnum:]]", "", names(dataset))
	row.names(dataset) = names(dataset)	
	dataset = round(dataset, digits=9)
	if (!openmx){
			require(sem, quietly=TRUE)
			a=b
			p = ""
			
			nms = unique(c(b[,1], b[,2]))
	
			#### write path equations
			for (i in 1:nrow(a)){
				if (a$Arrows[i] == 2) {symb="<->"} else {symb="->"}
				p = paste(p, a$From[i], symb, a$To[i], ", par", i, "\n", sep="")
			}
			end = which(nms %in% a$To)
			ex = which(!(nms %in% a$To))
			
			##### write endogenous variance equations
			if (length(nms)>0){
				for (i in 1:length(end)){
				p = paste(p, nms[end[i]], "<->", nms[end[i]], ", var[", nms[end[i]], "]\n", sep="")
				}
			}

			##### write exogenous variance equations			
			if (length(ex)>0){
				for (i in 1:length(ex)){
				p = paste(p, nms[ex[i]], "<->", nms[ex[i]], ", NA, 1\n", sep="")
				}
			}
			
			###### create folder to write to 
			cat(p, file="deleteme.txt")
			mod = specifyModel(file="deleteme.txt", quiet=TRUE)
			
			###### delete that folder
			cmd2 = paste("rm -rf ", getwd(), "/deleteme.txt", sep="")
			system(cmd2)
			se.m = suppressWarnings(sem(mod, S=as.matrix(dataset), N=1000, maxiter=10000, warn=FALSE))
			#attr(se.m, "class") = "DE.fit"
			se.m
	} else {
		require(OpenMx, quietly=TRUE)

			##### create connections between variables
			
			##### remove characters that give openmx trouble
		paths = mxPath(from=b[,1], to=b[,2], arrows=b[,3], values=b[,4], labels=paste(b[,1],"to", b[,2], "_", sep=""), free=T) 
	
			#### specify data as mx dataset
		#ord = order(row.names(dataset))		
		#dataset = dataset[ord, ord]
		data.set = mxData(dataset, type="cov", numObs=N)

			#### extract variable names
		variable.names = sort(unique(c(paths@from, paths@to)))
	
			##### find exogenous variables
		ex = which(!(variable.names %in% b$To))
	
			#### set all variances (if specified)
		if (fix.variances){	
			ar = rep(TRUE, times=length(variable.names))
			if (length(ex)>0) ar[ex] = FALSE	
			variances = mxPath(from=variable.names, arrows=2, connect="single", free=ar, values=1)
		} else {
			variances = mxPath(from=variable.names, arrows=2, connect="single", free=T, values=1)
		}
	
			#### put into one mxModel
		mxMod = mxModel("Randomly Generated Model", data.set, type="RAM", manifestVars=variable.names, variances, paths)
		
			#### don't estimate standard errors
		if (!calc.se){
			mxMod <- mxOption(mxMod, 'Calculate Hessian', 'No')
			mxMod <- mxOption(mxMod, 'Standard Errors', 'No')
		}
			#### run the model
		mxFit = mxRun(model = mxMod, silent=TRUE, suppressWarnings=TRUE)
		#attr(mxFit, "class") = "DE.fit"
		mxFit
	}
}
