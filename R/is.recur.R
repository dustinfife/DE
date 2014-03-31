is.recur = function(r.model, var.names){
	two.arrows = which(r.model$Arrows==2)
	if (length(two.arrows)>0){
		reformat.model = r.model[-two.arrows,]
	} else {
		reformat.model = r.model
	}
	reformat.model$From = as.character(reformat.model$From)
	reformat.model$To = as.character(reformat.model$To)
	for (o in 1:length(var.names)){
		froms = which(reformat.model$From == var.names[o])
		tos = which(reformat.model$To == var.names[o])
		reformat.model$From[froms] = o
		reformat.model$To[tos] = o
	}
	check = as.numeric(c(rbind(reformat.model$From, reformat.model$To)))
	recur = is.dag(graph(check))
	return(recur)	
}