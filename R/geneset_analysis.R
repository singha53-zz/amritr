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
#'library(dplyr)
#'data(pathways)

#'# gene list supplied by the user
#'set.seed(1)
#'genelist <- pathways %>%
#'  filter(Database == "KEGG") %>%
#'  dplyr::select(Genes) %>%
#'  unlist %>%
#'  as.character %>%
#'  sample(., 50)
#'
#'set.seed(2)
#'universe <- pathways %>%
#'  filter(Database == "KEGG") %>%
#'  dplyr::select(Genes) %>%
#'  unlist %>%
#'  as.character %>%
#'  sample(., 50) %>%
#'  c(., genelist) %>%
#'  unique
#'
#'# Perform gene set enrichment analysis
#'result <- gene_set_analysis(genes = genelist,
#'  genesets = pathways,
#'  study_universe = universe,
#'  min.set.size = 1)
#'
#'result$toptable %>% group_by(Type, Database) %>% slice(1:5)
#'
#' @export
gene_set_analysis = function(genes, genesets, study_universe, min.set.size = 1){
  ## Proportion of mapped background set
  prop.mapped.backgroundSet <- genesets %>%
    dplyr::group_by(Database) %>%
    summarise(prop.mapped.backgroundSet = round(100*length(intersect(unique(unlist(Genes)), study_universe))/length(study_universe), 1))

  ## unfiltered results
  genesets$universe <- genesets$Genes %>% unlist %>% unique %>% length
  unfiltered <- genesets %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(n_geneset = length(Genes)) %>%
    dplyr::filter(n_geneset > min.set.size) %>%
    dplyr::mutate(n_genes = length(genes),
      overlap = length(intersect(genes, Genes)),
      p_value = phyper(overlap - 1, n_genes,
        universe - n_genes, n_geneset, lower.tail = F)) %>%
    ungroup() %>%
    dplyr::group_by(Database) %>%
    dplyr::mutate(bhfdr = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(Database, p_value)
  unfiltered$Type <- "unfiltered"

  mapped.genes.unfiltered <- genesets %>%
    dplyr::group_by(Database) %>%
    summarise(mapped.genes.unfiltered = round(100*length(intersect(unique(unlist(Genes)), genes))/length(genes), 1))

  # filter genes using study_universe
  genesets <- genesets %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(Genes = list(intersect(Genes, study_universe)),
      n_geneset = length(Genes)) %>%
    dplyr::filter(n_geneset > min.set.size) %>% ungroup()
  genesets$universe <- genesets$Genes %>% unlist %>% unique %>% length

  filtered <- genesets %>%
    dplyr::mutate(n_genes = length(genes),
      overlap = length(intersect(genes, Genes)),
      p_value = phyper(overlap - 1, n_genes,
        universe - n_genes, n_geneset, lower.tail = F)) %>%
    ungroup() %>%
    dplyr::group_by(Database) %>%
    dplyr::mutate(bhfdr = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(Database, p_value)
  filtered$Type <- "filtered"
  toptable <- rbind(unfiltered, filtered)

  mapped.genes.filtered <- genesets %>%
    dplyr::group_by(Database) %>%
    summarise(mapped.genes.filtered = round(100*length(intersect(unique(unlist(Genes)), genes))/length(genes), 1))

  ## Summarise proportion of overlap
  prop <- cbind(mapped.genes.unfiltered, mapped.genes.filtered[,-1], prop.mapped.backgroundSet[,-1])

  return(list(toptable=toptable, prop=prop))
}


