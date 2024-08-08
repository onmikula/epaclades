#' Rooting `jplace` object.
#' 
#' @description
#' Roots the tree & adjusts edge labelling accordingly.
#'
#' @param jp an object of class `jplace`.
#' @param outgroup a character vector listing outgroup sequences or their unique identifier(s),
#'   e.g., a character string "outgroup", which is the default.
#' @returns A modified `jplace` object, i.e., its `tree` is rooted, `outgroup` component added to it,
#'   `placements$edge_num` is modified and if any placements were to the root edge,
#'   their rows in `placements` are modified too.
#' @export

root_jplace <- function(jp, outgroup="outgroup") {

	place <- jp$placements
	unrooted <- jp$tree
	outgroups <- unlist(lapply(outgroup, grep, x=unrooted$tip.label, value=TRUE))
	outgroups <- sort(unique(c(outgroups, unlist(lapply(outgroup, grep, x=unrooted$tip.label, value=TRUE, fixed=TRUE)))))
	rooted <- ape::root(unrooted, outgroup=outgroups, resolve.root=TRUE)	
	rooted <- ape::rotateConstr(rooted, rooted$tip.label)
	rooted$outgroup <- outgroups

	ubipart <- getbipartcodes(unrooted, rooted$tip.label)
	rbipart <- getbipartcodes(rooted, rooted$tip.label)
	rooted$edge.label <- unrooted$edge.label[match(rbipart, ubipart)]
	rootlabel <- rooted$edge.label[duplicated(rooted$edge.label)]
	place$edge_num <- rooted$edge.label[match(ubipart[match(place$edge_num, unrooted$edge.label)], rbipart)]
	rootplace <- place$edge_num == rootlabel
	if (any(rootplace)) {
		rplace <- place[rep(which(rootplace), each=2),]
		rplace$edge_num <- paste0(rplace$edge_num, c("L","R"))
		rplace$likelihood <- rplace$likelihood / 2
		rplace$like_weight_ratio <- rplace$like_weight_ratio / 2
		rplace$distal_length <- rplace$pendant_length <- NA
		if (!is.null(rplace$clade)) {
			rplace$clade <- NA
			rplace$position <- "ancestral"
		}
		place <- rbind(place[!rootplace,], rplace)
		place <- place[order(rownames(place)),]
	}
	basal <- rooted$edge.label == rootlabel
	rooted$edge.label[basal] <- paste0(rooted$edge.label[basal], c("L","R")[order(rooted$edge[basal,2])])

	jp$tree <- rooted
	jp$placements <- place	
	return(jp)

}



#' Bipartition codes.
#' 
#' @description
#' Assign codes to branches using the corresponding complementary sets of tips (bipartitions).
#'
#' @param phy an object of class `phylo`.
#' @param ref a character vector with tip labels (may or may not be tip labels of `phy`).
#' @returns A character vector with unique codes of tree branches, which are based on bipartitions defined by the branches.
#' @details The function is intended primarily as an internal function of `epaclades` package,
#'   so its input and output is tailored to the needs of [root_jplace].
#' @export

getbipartcodes <- function(phy, ref) {
	clades <- ape::prop.part(phy)
	labels <- attr(clades, "labels")
	clades <- cbind(diag(length(labels)), sapply(clades, function(x) as.numeric(ref %in% labels[x])))
	clades2 <- apply(clades[,phy$edge[,2]], 2, paste, collapse="")
	clades1 <- apply(abs(clades[,phy$edge[,2]]-1), 2, paste, collapse="")
	bipart <- apply(apply(cbind(clades1, clades2), 1, sort), 2, paste, collapse="-")
	return(bipart)
}
