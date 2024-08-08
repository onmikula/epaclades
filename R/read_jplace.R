#' Reading .jplace file.
#' 
#' @description
#' Reads an output of RAxML's EPA or pplacer in '.jplace' format.
#'
#' @param jp character, the name of '.jplace' output of EPA or pplacer.
#' @returns An object of class `jplace`, which is a list with components:
#'   `tree`, a `phylo` object with the unrooted reference tree used by EPA;
#'   `placements`, a data.frame with EPA output, i.e., list of possible placements
#'   of query sequences into the tree and their supports (likelihoods or posterior probabilities);
#'   and `software`, a character string indicating the software used (either `"pplacer"` or `"raxml"`).
#' @export

read_jplace <- function(jp) {

# input
	jp <- gsub("\\s", "", readLines(jp))

# tree
	tree <- paste0("(", gsub("^.*\"\\(|;.*$", "", jp[grep("\"tree\"", jp)]), ";")
	edge <- gregexpr("\\{[[:digit:]]+}", tree)
	edgelabel <- Map(paste0, regmatches(tree, edge), ":")
	regmatches(tree, edge) <- ""
	node <- gregexpr(":", tree)
	regmatches(tree, node) <- edgelabel
	tree <- ape::read.tree(text=tree)
	tiplabel <- regmatches(tree$tip.label, gregexpr("\\{\\d+}", tree$tip.label))
	tiplabel <- gsub("[\\{}]", "", unlist(tiplabel))
	nodelabel <- regmatches(tree$node.label, gregexpr("\\{\\d+}", tree$node.label))
	nodelabel <- lapply(nodelabel, function(x) gsub("[\\{}]", "", x))
	nodelabel <- unlist(ifelse(sapply(nodelabel, length) == 0, NA, nodelabel))
	edgelabel <- c(tiplabel, nodelabel)
	tree$tip.label <- gsub("\\{[[:digit:]]+}", "", tree$tip.label)
	tree$node.label <- NULL
	tree$edge.label <- edgelabel[tree$edge[,2]]

# placements
	place <- jp[grep("\"p\":", jp)]
	place <- setNames(gsub("^\\[|]$", "", regmatches(place, regexpr("\\[\\[.*]]", place))),
		gsub("\\[\"|\"]", "", regmatches(place, regexpr("\\[\".*\"]", place))))
	place <- lapply(place, function(p) unlist(regmatches(p, gregexpr("\\[.*]{1}?", p))))
	place <- lapply(place, function(p) gsub("^\\[|]$", "", unlist(p)))
	place <- lapply(place, function(p) do.call(rbind, lapply(strsplit(p, ","), as.numeric)))
	fields <- jp[grep("fields", jp):length(jp)]
	fields <- paste(fields[Reduce(":", grep("\\[|]", fields))], collapse="")
	fields <- unlist(strsplit(gsub("^.*\\[|].*$|\"", "", fields), ","))
	place <- setNames(data.frame(rep(names(place), sapply(place, nrow)), do.call(rbind, place)), c("id", fields))
	place$edge_num <- as.character(place$edge_num)

# software
	getstring <- function(pattern, text) regmatches(text, regexpr(pattern, text)) 
	soft <- jp[grep("\"invocation\"", jp)]
	soft <- unlist(lapply(c("raxml","pplacer"), getstring, text=tolower(soft)))

# output
	result <- list(tree=tree, placements=place, software=soft)
#	attr(result, "classified") <- FALSE
	class(result) <- "jplace"

	return(result)

}
