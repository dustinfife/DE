#'Outputs a RAM model as a graphic
#'
#'Outputs a pdf graphic of the model using graphviz
#'
#'@param model The RAM model to be imaged.
#'@param file location where the pdf file is to be stored
#'@param dotFilename not used.
#'@param rank.direction a string of the form, either "LR", or "TB", indicating
#'Left to right or top to bottom
#'@param max.rank a string of the form "V1, V2, V3", where variables 1-3 are
#'either at the most-right of the page or top-most (depending on whether
#'rank.direction = LR or TB)
#'@param min.rank same as above but will be in most-left or bottom-most
#'@param same.rank what variables should be lined up. Uses same format as
#'max.rank and min.rank.
#'@param delete.dot logical. Indicates whether the dot files should be deleted
#'(and thus only return pdf files)
#'@param silent Should R output information from system commands?
#'@param label The label attached to the graphviz diagram. It will be displayed
#'at the bottom of the image.
#'@param node.labels Should the nodes be labeled?
#'@author Dustin Fife
#'@examples
#'
#'# generate random model
#'rand.model = DE(LETTERS[1:5], paths=6)
#'# output to graphic
#'DE.graphviz(rand.model, min.rank="A, B, C")
#'# now go look at the root directory for the pdf graphic
#'@export
DE.graphviz <-
function (model, file="deleteme", dotFilename = "graphic", rank.direction="LR", max.rank=NULL, min.rank=NULL, same.rank=NULL, delete.dot=FALSE, silent=TRUE, label="", node.labels=TRUE) {

	#### open file connection    
	dot.file <- paste(file, ".dot", sep = "")
	handle <- file(dot.file, "w")
	on.exit(close(handle))
	graph.file <- paste(file, ".pdf", sep = "")    
	
	##### start writing the header of the file
	cat(file = handle, paste(""))
    cat(file = handle, paste("digraph \"", deparse(substitute(model)), "\" {\n", sep = ""))
    cat(file = handle, paste("  rankdir=", rank.direction, ";\n", sep = ""))
    cat(file = handle, paste("  size=\"8,8\";\n", sep = ""))
    cat(file = handle, paste("  node [fontname=\"Times-Roman\" fontsize=8 shape=box];\n", sep = ""))
    cat(file = handle, paste("  node [fontname=\"Times-Roman\" fontsize=14];\n", sep = ""))    
    cat(file = handle, "  center=FALSE;\n")
	
	##### specify the rank
	if (!is.null(min.rank)) {
        min.rank <- paste("\"", min.rank, "\"", sep = "")
        min.rank <- gsub(",", "\" \"", gsub(" ", "", min.rank))
        cat(file = handle, paste("  {rank=min ", min.rank, "}\n", 
            sep = ""))
    }
    if (!is.null(max.rank)) {
        max.rank <- paste("\"", max.rank, "\"", sep = "")
        max.rank <- gsub(",", "\" \"", gsub(" ", "", max.rank))
        cat(file = handle, paste("  {rank=max ", max.rank, "}\n", 
            sep = ""))
    }
    if (!is.null(same.rank)) {
        for (s in 1:length(same.rank)) {
            same <- paste("\"", same.rank[s], "\"", sep = "")
            same <- gsub(",", "\" \"", gsub(" ", "", same))
            cat(file = handle, paste("  {rank=same ", same, "}\n", 
                sep = ""))
        }
    }
  

	###### start writing the paths
	if (node.labels){
		direction = model$Arrows
	    for (par in 1:nrow(model)) {
	        if ((direction[par] == 1)) { direction[par] = ""} else { direction[par] = " dir=both"}
	            cat(file = handle, paste("  \"", model$From[par], 
	                "\" -> \"", model$To[par], "\" [label=\"", 
	                round(model$Values[par], digits=3), "\"", direction[par], "];\n", sep = ""))
	    }
    } else {
		direction = model$Arrows
	    for (par in 1:nrow(model)) {
	        if ((direction[par] == 1)) { direction[par] = ""} else { direction[par] = " [dir=both]"}
	            cat(file = handle, paste("  \"", model$From[par], 
	                "\" -> \"", model$To[par], "\"", direction[par], ";\n", sep = ""))
	    }    	
    }

	##### put title
	cat(file = handle, paste("  lobelloc=\"t\";\nlabel=\"", label, "\";\n", sep=""))

    
    ##### output dot to a file, produce pdf output, then delete (optional) the dot file
    cat(file = handle, "}\n")
    if (!missing(file)) {
        cmd <- paste("dot -T", "pdf", " -o ", graph.file, 
            " ", dot.file, sep = "")
        if (!silent) cat("Running ", cmd, "\n")
        result <- try(system(cmd))
        if (delete.dot){
	        cmd2 = paste("rm ", dot.file, sep="")
	        if (!silent) cat("Running ", cmd2, "\n")
    	    try(system(cmd2))
    	}    
    }
}
