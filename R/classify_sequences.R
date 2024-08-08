#' Classification of sequences.
#' 
#' @description
#' Classifies sequences to clades (or bipartitions).
#'
#' @param jp an object of class `jplace` pre-treated by `classify_jplace`,
#'   i.e., containing classification to clades / bipartitions.
#' @param keep numeric, no. of placements per clade to keep, corresponds to 'epa-­keep-­placements' in RAxML. 
#'   `keep=1` means the maximum likelihood placement, while `keep=NULL` (default) means all placements.
#' @param probthr numeric, likelihood weight threshold, corresponds to 'epa­-prob­-threshold' in RAxML. 
#' @param cumprobthr numeric, cummulative likelihood weight threshold, corresponds to
#'    'epa­-accumulated­-threshold'in RAxML.
#' @returns A data frame with components:
#'   `id` with labels of query sequences;
#'   `clade` with labels of clade / bipartition;
#'   `position` which is either `"crown"` or `"stem"`;
#'   `probability` with probability of classification, i.e., sum of probabilities over the branches
#'   belonging to the clade / bipartition. In this context, 'probability' means either posterior probability
#'   or relative likelihood weight, i.e., a quantity which sums up to unity over all branches of the tree.
#'   `n_branch` which is the number of branches which contributed to the sum.
#' @export

classify_sequences <- function(jp, keep=NULL, probthr=0, cumprobthr=1) {
	if (is.null(jp$placements$clade)) {
		stop("the 'jplace' object was not pre-treated by 'classify_jplace'")
	}
	if (is.null(keep)) {
		keep <- ape::Nedge(jp$tree)
	}
	probability <- c(pplacer="post_prob", raxml="like_weight_ratio")[jp$software]
	place <- jp$placements
	place$label <- paste(place$clade, place$position, sep="_")
	place <- split(place, place$id)
	for (i in seq_along(place)) {
		iprob <- place[[i]][,probability]
		ikeep <- which(iprob >= probthr & cumsum(iprob) <= cumprobthr)
		if (length(ikeep) > 0) {
			ikeep <- seq(min(c(max(ikeep), keep)))
		} else {
			ikeep <- FALSE
		}
		place[[i]] <- place[[i]][ikeep,,drop=FALSE]
		if (nrow(place[[i]]) > 0) {
			place[[i]] <- split(place[[i]], place[[i]]$label)
			for (j in seq_along(place[[i]])) {
				prob <- data.frame(probability=sum(place[[i]][[j]][,probability]), n_branch=nrow(place[[i]][[j]]))
				place[[i]][[j]] <- cbind(place[[i]][[j]][1, c("id","clade","position"), drop=FALSE], prob)
			}
			place[[i]] <- do.call(rbind, place[[i]])
			place[[i]] <- place[[i]][order(place[[i]]$probability, decreasing=TRUE),,drop=FALSE]
		}
	}
	place <- do.call(rbind, place)
	rownames(place) <- NULL
	return(place)
}



#' Best supported placement.
#' 
#' @description
#' Finds the maximum likelihood or maximum aposteriori placement of each sequence.
#'
#' @param place A data frame, the output of `classify_sequences`.
#' @returns A reduced data frame with best supported classification for each query sequence.
#'   (ML classification for EPA, MAP for pplacer.)
#' @export
best_placement <- function(place) {
	bestplace <- place[match(unique(place$id), place$id),]
	bestplace$totalprob <- bestplace$probability
	for (i in seq(nrow(bestplace))) {
		matchid <- place$id == bestplace$id[i]
		matchclade <- place$clade == bestplace$clade[i]
		bestplace$totalprob[i] <- sum(place$probability[which(matchid & matchclade)])
	}
	return(bestplace)
}
