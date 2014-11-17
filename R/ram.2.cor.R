##' Generate a correlation matrix from a Reticular Action Model (RAM)
##'
##' Generate a correlation matrix from a RAM
##'	
##' Although generating an implied covariance matrix from a RAM model is easy, doing so
##' for a correlation is not as trivial since the diagonals are constrained to one. The symmetric matrix (S)
##' must be formulated in such as way that the sum of the variances equals one. Although there may exist a
##' closed-form solution, I know of none. So, instead, I simply used the optim function to do it. The user 
##' specifies the path coefficients in the form of a RAM matrix and the algorithm estimates what the residual
##' correlations ought to be.
##' @param RAM a Reticular Action model with columns "From", "To", "Arrows", and "Values"
##' @references An openmx forum about this: http://openmx.psyc.virginia.edu/thread/3866
##' @return A correlation matrix generated from a RAM model
##' @author Dustin Fife
##' @export
##' @examples
##' ## generate a simple exercise model
##' RAM = data.frame(matrix(c(
##' 	"Age", "Energy", 1, .1,
##' 	"Age", "Self-Efficacy", 2, .1,
##' 	"Energy", "Exercise", 1, .6,
##' 	"Eating", "Energy", 1, .4,
##' 	"Eating", "Self-Efficacy", 1, .3,
##' 	"Age", "Eating", 2, .15,
##' 	"Self-Efficacy", "Exercise", 1, .5), ncol=4, byrow=T))
##' names(RAM) = c("From", "To", "Arrows", "Values")	
##' ram.2.cor(RAM)
ram.2.cor = function(RAM){
	
	#### first create function that can be optimized
	opfunc = function(diag=c(.5, .5, .5, .5), RAM, op=TRUE){
		
		#### extract the names of all the variables
		all.vars = as.character(unique(unlist(c(RAM[,1:2]))))
		
		#### create asymmetric matrix
		mA = data.frame(matrix(0,nrow=length(all.vars), ncol=length(all.vars)), row.names=all.vars)
		names(mA) = c(all.vars)
		for (j in 1:nrow(RAM)){
			col = which(names(mA) == RAM$From[j] & RAM$Arrows[j] == 1)
			rw = which(names(mA) == RAM$To[j] & RAM$Arrows[j] == 1)		
			mA[rw, col] = as.numeric(as.character(RAM$Values[j]))
		}
	
		##### create symmetric matrix (temporarily fill in endogenous variances)
		end = names(mA)[which(rowSums(mA)>0)]
		mS = data.frame(diag(x=diag, length(all.vars)), row.names=all.vars)
		names(mS) = row.names(mS)
		for (i in 1:nrow(RAM)){
			if (RAM$Arrows[i] == 2){
				col = which(names(mA) == RAM$To[i])
				row = which(names(mA) == RAM$From[i])
				mS[col,row] = as.numeric(as.character(RAM$Values[i]))
				mS[row,col] = as.numeric(as.character(RAM$Values[i]))
			}
		}		
		
		#### create factor and identity matrices
		mI = diag(1, nrow(mA))
		mF = mI
		
		#### compute implied covariance matrix
		C = (solve(mI-mA) %*% data.matrix(mS) %*% t(solve(mI-mA)))

		#### return the sum of squares of the diagonals
		diags = sum((diag(C) - 1)^2)
		
		if (op){
			return(diags)
		} else {
			return(C)
		}
	}
	
	### give a note
	cat("Doing a brute-force search for the correct correlations. Give me a minute.\n\n")
	
	#### now we optimize it
	vars = as.character(unique(unlist(c(RAM[,1:2]))))
	op.res = optim(rep(.5, times=length(vars)), opfunc, RAM=RAM)
	
	#### and return the original matrix
	RAM.returned = round(opfunc(diag=op.res$par, RAM, op=FALSE), digits=2)
	RAM.returned
	
}
