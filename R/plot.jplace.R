#' Plotting phylogenetic placements.
#' 
#' @description
#' A plot method for class `jplace`, it highlights branches with placements of query sequences.
#'
#' @param x an object of class `jplace`.
#' @param subset character vector with sequence names.
#' @param keep numeric, how many best placements to keep for display.
#' @param prob numeric, ummulative probability threshold up to which the placements are displayed.
#' @param type character string indicating the type of tree, reasonable is `"radial"` (default) and `"unrooted"`.
#' @param col color to highlight branches with placements.
#' @param bg color of background branches.
#' @param file optional character string, the name of .pdf file into which the plot is saved.
#' @param ... other parameters of [ape::plot.phylo].
#' @export

plot.jplace <- function(x, subset=NULL, keep, prob=1.00, type="radial", col=3, bg=1, file, ...) {
	probability <- c(pplacer="post_prob", raxml="like_weight_ratio")[x$software]
	place <- x$placements
	if (!probability %in% names(place)) {
		stop("edges of the tree in 'jplace' object are yet to be classified into clades")
	}
	if (!is.null(subset)) {
		place <- place[place$id %in% subset,,drop=FALSE] 		
	}
	place <- split(place, place$id)
	if (prob < 1) {
		place <- lapply(place, function(p) p[1:min(which(cumsum(p[,probability]) >= prob)),,drop=FALSE])
	}
	if (!missing(keep)) {
		place <- lapply(place, function(p) p[1:min(keep,nrow(p)),,drop=FALSE])		
	}
	place <- do.call(rbind, place)
	edgecol <- c(bg, col)[x$tree$edge.label %in% place$edge_num + 1]
	if (!missing(file)) {
		pdf(file)
	}	
	plot(x$tree, type=type, edge.color=edgecol, ...)
	if (!missing(file)) {
		invisible(dev.off())
	}	
}


#' Branch labels.
#' 
#' @description
#' Makes labels indicating which query sequences are most probably placed on given branches.
#'
#' @param jp An object of class `jplace`.
#' @param classification A data frame, output of `classify_sequences`.
#' @param sep A character string separating labels of query sequences placed on the same branch.
#' @returns A character vector with concatenated id's of query sequences most probably placed on given branch.
#'   The vector elements are named by branch labels from `x$tree$edge.label`.
#' @export
make_plotlabels <- function(jp, classification, sep="\n") {
	best <- classification[match(unique(classification$id), classification$id),]
	probability <- c(pplacer="post_prob", raxml="like_weight_ratio")[jp$software]
	place <- jp$placements
	place <- split(place, place$id)
	best <- best[match(names(place),best$id),]
	for (i in seq_along(place)) {
		place[[i]] <- place[[i]][place[[i]]$clade == best$clade[i] & place[[i]]$position == best$position[i],,drop=FALSE]
		place[[i]] <- place[[i]][which.max(place[[i]][,probability]),,drop=FALSE]
	}
	place <- do.call(rbind, place)
	edges <- split(place$id, place$edge_num)
	edgelab <- sapply(edges, paste, collapse=sep)
	edgelab <- setNames(edgelab, match(names(edgelab), jp$tree$edge.label))
	return(edgelab)
}

