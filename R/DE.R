#'Generates a random model using RAM specification.
#'
#'This function generates a single random path model.
#'
#'Making restrictions is simple. For example, suppose one variable is Age.
#'Obviously, Age should probably not be endogenous, so the user can specify Age
#'as an endogenous variable. That is done by creating a matrix where the
#'columns correspond to "From", "To", and "Include." For example, to specify
#'that A must cause B, one would insert in the first row of the matrix c("A",
#'"B", "1"). To specify that nothing can cause a variable (i.e., to make a
#'variable exogenous), one would leave the "From" column as "". For example,
#'the Age example would have c("", "Age", "0").
#'
#'@param variable.names A vector of variable names.
#'@param paths The number of paths to be randomly generated. Can be a single value or a vector of two integers specifying a range. 
#'@param restrictions What kind of restrictions are set. (See details).
#'@param prop.arrows What proportion of exogenous variables should be correlated. Defaults to .2.
#'@param allow.orphaned Should orphaned variables be allowed when random models
#'are generated?
#'@param allow.bidir Should bidirectional arrows be allowed? (Note: this is not
#'the same as a correlation).
#'@param corr.exogenous Should exogenous variables be correlated?
#'@return Returns a RAM matrix.
#'@author Dustin Fife
#'@seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#'@references Fife, D.A., Rodgers, J.L., & Mendoza, J. L. (2013). Model
#'conditioned data elasticity in path analysis: Assessing the "confoundability"
#'@export
#'@examples
#'
#'rest = matrix(c("A", "B", "1",
#'						"", "A", "0"), nrow=2, byrow=TRUE)
#'DE(variable.names=LETTERS[1:6], paths=c(6,7), restrictions = NULL, prop.arrows = 0.2, allow.orphaned=FALSE, allow.bidir=FALSE, allow.cov.endogenous=FALSE)
#'


####
#variable.names=LETTERS[1:6]; paths=8; restrictions=NULL; prop.arrows=.3; allow.orphaned=F; allow.bidir=F; allow.cov.endogenous=F; corr.exogenous=T
DE <-
function(variable.names, paths, restrictions=NULL, prop.arrows=.2, allow.orphaned=FALSE, allow.bidir=FALSE, corr.exogenous=FALSE){
	
	variables = length(variable.names)
	
	#### alerts
	if (length(paths)==1){
		if (paths<variables){
			stop("The number of paths must be at least as great as the number of variables.")
		}
	} else {
		cond = (paths<variables)
		if (cond[1] | cond[2]){
			stop("The number of paths must be at least as great as the number of variables.")
		}
	}
	if (length(variable.names) != variables){
		stop(paste("The object 'variable.names' must have exactly ", variables, " elements.", sep=""))
	}
	
	#### randomly sample the number of paths, given constraints
	if (length(paths)>1){
		paths = round(runif(1, paths[1]-.5, paths[2]+.5))
	}

	
	#### come up with all possible ways of connecting variables
	poss = expand.grid(variable.names, variable.names) 
	poss$must = 0 
	
	#### remove variances
	poss = poss[-which(poss[,1] == poss[,2]),] 

	#### check constraints
	if (!is.null(restrictions)){
		if (ncol(restrictions) != 3){
		stop("Your 'restrictions' matrix must have three columns: From, To, and Include.")
	}}
	

	#### remove/add user-specified constraints
	if (!is.null(restrictions)){
		for (i in 1:nrow(restrictions)){  
			
			#### take care of whether a variable is specified as endogenous or not
			if (restrictions[i,1] == "" | restrictions[i,2] == ""){
				
				#### make sure they don't have both empty
				if (restrictions[i,1] == "" & restrictions[i,2] == "") {
					stop("You have a row where more than one slot is empty.")
				}
				
				#### identify the variables and where they correspond
				column = which(restrictions[i,1:2] != "")
				#### find rows to delete (to avoid mutual causation)
				column2 = which(restrictions[i,1:2] == "")				
				rws2 = which(poss[,column] == restrictions[i,column])
				delRows = which(poss[,column2] == restrictions[i,column])

				##### include them
				if (restrictions[i,3] == 1){
					poss[rws2,3] = 1
					### remove bidirectional
					if (length(delRows)>0) poss = poss[-delRows,]
				} else {  ### or remove them
					poss = poss[-rws2,]
				}
				##### take care of other restrictions
			} else {
				rws = which(as.character(poss[,1]) == restrictions[i,1] & 
							as.character(poss[,2]) == restrictions[i,2])
				rwsDel = which(as.character(poss[,2]) == restrictions[i,1] & 
							as.character(poss[,1]) == restrictions[i,2]) 			
						
				if (restrictions[i,3] == 1){
					poss[rws,3] = 1
					### remove bidirectional
					if (length(rwsDel)>0) poss = poss[-rwsDel,] 
				} else {
					poss = poss[-rws,]
				}
			}
		}
	}

	
	###### TODO: tell user how many possible models there are

	###### first make sure there's no "orphaned" variables
	if (!allow.orphaned){
	random.model = matrix(nrow=paths, ncol=2)
	for (i in 1:variables){
		
		#### see if that variable has a one next to it
		subset.poss = subset(poss, Var1==variable.names[i] | Var2==variable.names[i])
		#### if the variable doesn't yet have a one

		if (!(1 %in% subset.poss$must)){
			
			##### randomly select of the remaining
			remaining = which(poss$Var1==variable.names[i] | poss$Var2 == variable.names[i])
			rand.samp = sample(remaining, 1)
			
			##### make it a one
			poss$must[rand.samp] = 1
			##### remove bidirectional
			r = which(poss$Var2 == poss[rand.samp,"Var1"] & poss$Var1 == poss[rand.samp,"Var2"])
			if (length(r)>0) poss = poss[-r,]			
		}
			#### randomly select one of those paths
	}}



	##### for the non-orphaned variables, randomly select remaining paths
	used = length(which(poss$must == 1))
	left = paths-used
	if (left<0){
		stop("The user-specified constraints cannot be met with the number of paths 
				you have selected. Either decrease the number of constraints or increase 
				the number of allowable paths.")
	} else if (left>0){
		for (i in 1:left) {
			###### randomly sample "left" paths from remaining, eliminating bidirectional			
			rem.rws = which(poss$must==0)
			rand.row = sample(rem.rws, 1)
			#### select that row
			poss[rand.row,3] = 1
			#### remove its twin (if specified)
			if (!allow.bidir){
				del.row = which(poss[,2] == poss[rand.row,1] & poss[,1] == poss[rand.row,2])
				if (length(del.row)>0) poss = poss[-del.row,]
			}
		}
	}
	random.model = poss[which(poss$must==1),1:2]


	##### determine number of arrows
	ar = rep(NA, times=nrow(random.model))

	#### fix user specified arrows
	if (!is.null(restrictions)){	
	for (i in 1:nrow(restrictions)){
		rw = which(random.model[,1] == restrictions[i,1] & random.model[,2] == restrictions[i,2])
		ar[rw] = as.numeric(restrictions[i,3])
	}}
	

	
	#### randomly specify remaining arrows
	# ar[is.na(ar)] = runif(length(which(is.na(ar))))
	# ar[ar<prop.arrows] = 2
	# ar[!(ar==1 | ar==2)] = 1
	ar[is.na(ar)] = 1

	##### set "values" column randomly
	vals = .2#runif(paths, 0, .4)
	
	##### combine into mxReady object
	random.model = data.frame(random.model)
	names(random.model) = c("From", "To")
	random.model$Arrows = ar
	random.model$Values = vals
	random.model$From = as.character(random.model$From)
	random.model$To = as.character(random.model$To)	

	# ##### eliminate the possibility of an endogenous covariance	
	# if (!allow.cov.endogenous){
		# ####### determine if endogenous variable has bidirectional arrow
		# #### find those that are endogenous for sure
		# sb = subset(random.model, Arrows==1)
		
		# #### if any row has a single arrow
		# if (nrow(sb)<nrow(random.model)){
			# ##### loop through every row 
			# for (k in 1:nrow(random.model)){
				# #### select two-headed arrows
				# sb2 = subset(random.model, Arrows==2)				
				# #### select one-headed arrows
				# sb = subset(random.model, Arrows==1)		
						
				# #### look at those who have single arrows going into them (endogenous)
				# end = variable.names[which(variable.names %in% sb$To)]				
				# s.end = end[which(!(end %in% sb$From))]
				
				# #### if the current row is double-headed
				# if (random.model$Arrows[k]==2){
					# ### see if either the from or the to row is endogenous
					# if (sum(end %in% random.model$To[k] | end %in% random.model$From[k])>0){
						
						# # #### if both are strinctly endogenous, leave it as 2
						# # if (sum(s.end %in% random.model$To[k])>0 & sum(s.end %in% random.model$From[k])>0){
							# # random.model$Arrows[k] = 2
						# # } else {
						# #### if one is endogenous and the other exogenous, set to one
							# random.model$Arrows[k] = 1
						# #}
					# }
				# }	#### if arrows is two
			# }	### loop through all rows
		# }		### condition to correct endogeneity
	# }
	
	##### make exogenous variables correlated (if specified)
	if (corr.exogenous){
		sb = subset(random.model, Arrows==1)
		
		##### find endogenous/exogenous variables
		end = variable.names[which(variable.names %in% sb$To)]				
		s.end = variable.names[which(!(variable.names %in% end))]
		
		##### if there's more than one exogenous variable
		if (length(s.end)>1){
			
			### come up with all possible ways of connecting endogenous variables
			allofem = data.frame(t(combn(s.end, 2))); names(allofem) = c("From", "To")
			
			### randomly select which rows have arrows
			randnum = runif(nrow(allofem))
			whichrows = randnum<prop.arrows
			
			msg = paste(length(s.end), "exogenous variables with", sum(whichrows), "being double-headed.\n")
			#print(msg)
				
			### fill out rest of allofem
			allofem$Arrows=1
			allofem$Values=.2

			#### fill in double headed arrows
			if (sum(whichrows)>0){
				allofem$Arrows[whichrows] = 2
			}
			
			#print(allofem[allofem$Arrows==2,])
			
			#### remove those that are not double-headed
			#allofem = allofem[-which(!(allofem$Arrows==2)),]

			#### put back with rest of model			
			random.model = rbind(random.model, allofem[allofem$Arrows==2,])
		}
		
	}
	
	###### return model
	return(random.model)
}

s.end = LETTERS[1:3]