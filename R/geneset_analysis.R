#' @title Gene Set Analysis (GSA)
#' @description
#' \code{gene_set_analysis} Perform gene set analysis using the hypergeometric test
#' @param genelist vector of gene symbols
#' @param genesets list of gene sets
#' @param universe background list of genes: all the genes used in the study
#' @details
#' This function uses the hypergeometric function to test for over-representation of genesets in a user-specified list of genes. The unfiltered genesets contain all genes that are part of certain databases, whereas filtered genesets only contain genes that were used in the experiment (in order prevent sampling bias, \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0761-7}).
#' @return The output is a dataframe with 4 columns (gene sets (genesets), p-value (pval), bhfdr (BH-FDR), Type (unfiltered or filtered)
#' @examples
#'
#' data(pathways)
#'
#' # genelist and KEGG pathways (gene sets)
#' set.seed(1)
#' genelist <- sample(unlist(unlist(kegg)), 50)
#' set.seed(2)
#' universe <- unique(c(genelist, sample(unlist(unlist(kegg)), 300)))
#'
#' # Perform gene set enrichment analysis
#' result <- gene_set_analysis(genelist, genesets, universe)
#'
#' result %>% group_by(Type) %>% slice(1:5)
#' @export
gene_set_analysis = function(genelist,
                             genesets,
                             universe){
# unfiltered
unfiltered <- hypertest(genelist, genesets)
unfiltered$Type <- "unfiltered"

## filtered
filtered_genesets <- lapply(genesets, function(i) intersect(i, universe))
filtered <- hypertest(genelist, filtered_genesets)
filtered$Type <- "filtered"

## print percentage of mapping and unmapped identifiers in the unfiltered genelist
cat(paste0("% of mapped identifiers in unfiltered genesets = ",
  100*round(length(intersect(genelist, unlist(genesets)))/length(genelist)), "%"), fill = TRUE)
cat(paste0("% of unmapped identifiers in unfiltered genesets = ",
  100*round(length(setdiff(genelist, unlist(genesets)))/length(genelist)), "%"), fill = TRUE)

## print percentage of overlapping and non-overlapping identifiers in the background set
cat(paste0("% of mapped identifiers in the background set = ",
  100*round(length(intersect(universe, unlist(genesets)))/length(universe)), "%"), fill = TRUE)
cat(paste0("% of unmapped identifiers in the background set = ",
  100*round(length(setdiff(universe, unlist(genesets)))/length(universe)), "%"), fill = TRUE)

rbind(unfiltered, filtered)
}


#' @title hypergeometric Test
#' @description
#' \code{hypertest} performs hypergeometric tests across a list of gene sets
#' @param genelist vector of gene symbols
#' @param genesets list of gene sets
#' @details
#' This function uses the hypergeometric function to test for over-representation of genesets in a user-specified list of genes.
#' @return The output is a dataframe with 3 columns (gene sets (genesets), p-value (pval), and bhfdr (BH-FDR)
#' @export
hypertest = function(genelist, genesets){
  data.frame(pval = sapply(genesets, function(i){
    uni <- unlist(genesets)
    input <- genelist
    n_input = length(input)
    n_geneset = length(i)
    intersect = length(intersect(input, i))
    stats::phyper(intersect - 1, n_input,
      length(uni) - n_input, n_geneset, lower.tail = F)
  })) %>%
    mutate(genesets = rownames(.),
      bhfdr = p.adjust(pval, "BH")) %>%
    arrange(pval) %>%
    dplyr::select(genesets, pval, bhfdr)
}
