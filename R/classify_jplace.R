#' Classification of `jplace` object.
#' 
#' @description
#' Classifies tree edges with any phylogenetic placements to clades (or bipartitions).
#'
#' @param jp an object of class `jplace`.
#' @param clades a matrix or data frame defining classification of tips (1st column) to clades / bipartitions
#'   (2nd column). The classification may not be comprehensive (not including all the tips), but the clades
#'   must not be overlapping or mutually nested.
#' @param ancestral Logical, whether to classify ancestral branches; The default option is `TRUE`,
#'   but it is relevant only if the tree is rooted or the classification to clades is comprehensive
#'   (otherwise it is switched to `FALSE` and issues a warning).
#' @returns A object of class `jplace` with additional components `clade` and `position`
#'   in the `placements` data frame, which indicate classification of the branches.
#' @export

classify_jplace <- function(jp, clades, ancestral=TRUE) {

	tree <- jp$tree
	tree <- ape::rotateConstr(tree, tree$tip.label)
	place <- jp$placements
	rooted <- ape::is.rooted(tree)
	comprehensive <- all(tree$tip.label %in% clades[,1])
	clades <- setNames(as.data.frame(clades[,1:2]), c("tip", "clade"))
	clades <- clades[clades$tip %in% tree$tip.label,]
	clades <- lapply(split(clades$tip, clades$clade), match, table=tree$tip.label)
	
	if (isTRUE(ancestral) & isFALSE(rooted | comprehensive)) {
		ancestral <- FALSE
		message <- c("the classification is not comprehensive", "the tree is not rooted")[c(!comprehensive, !rooted)]
		if (length(message) > 1) {
			message <- paste(message[1], "and", message[2])	
		}
		message <- paste0(toupper(substr(message, 1, 1)), substr(message, 2, nchar(message)), ", 'ancestral' argument was set to FALSE.")
		warning(message)
	}

	monophyletic <- sapply(clades, function(x) ape::is.monophyletic(tree, tips=x))
	if (all(!monophyletic)) {
		stop("none of the clades is monophyletic in the tree")
	} else if (any(!monophyletic)) {
		nonmono <- names(clades)[! monophyletic]
		clades <- clades[monophyletic]
		single <- length(nonmono) == 1
		if (!single) nonmono <- paste(paste(nonmono[-length(nonmono)], collapse=", "), "and", nonmono[length(nonmono)])
		warning(paste("The", c("clades", "clade")[single+1], nonmono, c("are", "is")[single+1], "not monophyletic in the tree."))
	}

	if (length(clades) > 1) {
		overlapping <- matrix(FALSE, length(clades), length(clades))
		for (i in 2:length(clades)) {
			for (j in 1:(i-1)) {
				overlapping[i,j] <- length(intersect(clades[[i]], clades[[j]])) > 0
			}
		}	
		if (any(overlapping)) {
			stop("some of the clades are overlapping")
		}
	}

	branches <- lapply(clades, bipartition, phy=tree)
	mrcas <- lapply(branches, function(b, phy) intersect(phy$edge[b[,"crown"],], phy$edge[b[,"stem"],]), phy=tree)
	crowns <- lapply(branches, function(b) tree$edge.label[b[,"crown"]])
	crowns <- data.frame(edge=unlist(crowns), clade=rep(names(crowns), sapply(crowns, length)), position="crown", row.names=NULL)	
	stems <- lapply(branches, function(b) tree$edge.label[b[,"stem"]])
	stems <- data.frame(edge=unlist(stems), clade=rep(names(stems), sapply(stems, length)), position="stem", row.names=NULL)
	duplstems <- stems$edge %in% stems$edge[duplicated(stems$edge)]
	if (any(duplstems)) {
		duplbranch <- tree$edge.label %in% stems$edge[duplstems]
		for (cl in unique(stems$clade[duplstems])) {
			branches[[cl]][duplbranch,"stem"] <- FALSE
		}
		stems <- stems[!duplstems,]
	}

	if (isTRUE(ancestral) & isFALSE(rooted)) {
		lastancestors <- function(b, a, phy) ifelse(any(b[,"stem"]), phy$edge[b[,"stem"],1], a)
		ancnodes <- mapply(lastancestors, branches, mrcas, MoreArgs=list(phy=tree), SIMPLIFY=FALSE)
		ancnodes <- sort(unlist(ifelse(sapply(ancnodes, length) == 0, clades, ancnodes)))
		edge <- apply(apply(tree$edge, 1, sort), 2, paste, collapse="-")
		ancedges <- rep(FALSE, length(edge))
		while (length(ancnodes) > 1) {
			path <- ape::nodepath(tree, ancnodes[1], ancnodes[length(ancnodes)])
			ancnodes <- c(ancnodes[1], setdiff(ancnodes, path))
			path <- cbind(path[-1], path[-length(path)])
			path <- apply(apply(path, 1, sort), 2, paste, collapse="-")
			ancedges <- ancedges | edge %in% path
		}
		ancestors <- data.frame(edge=tree$edge.label[ancedges], clade=NA, position="ancestral", row.names=NULL)
	} else if (isTRUE(ancestral)) {
		root <- ape::Ntip(tree) + 1
		bottom <- sapply(branches, function(b, phy) phy$edge[b[,"stem"],1], phy=tree)
		edge <- setNames(apply(apply(tree$edge, 1, sort), 2, paste, collapse="-"), tree$edge.label)
		ancedges <- lapply(bottom, function(to, phy, root) ape::nodepath(phy, from=root, to=to), phy=tree, root=root)
		ancedges <- lapply(ancedges, function(path) cbind(path[-1], path[-length(path)]))
		ancedges <- lapply(ancedges, function(path) apply(apply(path, 1, sort), 2, paste, collapse="-"))
		ancedges <- sapply(ancedges, function(path, edge) edge %in% path, edge=edge)
		rownames(ancedges) <- names(edge)
		ancedges <- ancedges[apply(ancedges, 1, any),,drop=FALSE]
		cladelabels <- apply(ancedges, 1, function(x, lab) paste0("anc_", paste(lab[x], collapse="-")), lab=names(clades))
		ancestors <- data.frame(edge=rownames(ancedges), clade=cladelabels, position="ancestral", row.names=NULL)
	} else {
		ancestors <- data.frame(edge=character(0), clade=character(0), position=character(0))
	}
	
	edges <- do.call(rbind, list(crowns, stems, ancestors))
	ii <- match(place$edge_num, edges$edge)
	place$clade <- edges$clade[ii]
	place$position <- edges$position[ii]
	jp$placements <- place
	return(jp)

}


#' Branches of a bipartition.
#' 
#' @description
#' Finds tree branches belonging to a bipartition defined by the specified collection of tips.
#'
#' @param tips a numeric vector with indices of tree tips belonging to the bipartition.
#' @param phy an object of class `phylo`.
#' @returns A matrix of mode `logical` whose columne correspond to tree branches and its two columns indicate,
#'   whether the branch belongs to the crown of the bipartition and the second whether it belongs to its stem.
#' @details The function is intended primarily as an internal function of `epaclades` package,
#'   so its input and output is tailored to the needs of [classify_jplace].
#' @export

bipartition <- function(tips, phy) {
	if (length(tips) > 1) {
		anc <- ape::getMRCA(phy, tips)
		paths <- lapply(tips, function(tip) ape::nodepath(phy, from=anc, to=tip))
		paths <- lapply(paths, function(p) cbind(p[-length(p)], p[-1]))
		paths <- lapply(paths, function(p) apply(apply(p, 1, sort), 2, paste, collapse="-"))
		edges <- apply(apply(phy$edge, 1, sort), 2, paste, collapse="-")
		crown <- edges %in% unlist(paths)
		incl <- as.numeric(phy$edge[crown,])
		stem <- xor(phy$edge[,1] %in% incl, phy$edge[,2] %in% incl)
		bipartmat <- cbind(crown=crown, stem=stem)
	} else {
		bipartmat <- cbind(crown=FALSE, stem=phy$edge[,2] %in% tips)
	}
	return(bipartmat)
}
