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
gene_set_analysis = function(genes, genesets, study_universe, alpha, topranked = 5, min.set.size = 1){
  ## unfiltered results
  genesets$universe <- genesets$Genes %>% unlist %>%
    unique %>% length
  unfiltered <- genesets %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(n_geneset = length(Genes)) %>%
    dplyr::filter(n_geneset > min.set.size) %>%
    dplyr::mutate(n_input = length(genes),
      overlap = length(intersect(genes, Genes)),
      p_value = phyper(overlap - 1, n_geneset,
        universe - n_geneset, n_input, lower.tail = F)) %>%
    ungroup() %>%
    dplyr::group_by(PathwayType) %>%
    dplyr::mutate(bhfdr = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(PathwayType, p_value)
  subset.unfiltered <- summarise(unfiltered, sig=sum(bhfdr < alpha)) %>%
    filter(sig != 0) %>% dplyr::select(PathwayType)
  unfiltered <- rbind(filter(unfiltered, PathwayType %in% subset.unfiltered$PathwayType) %>% filter(bhfdr < alpha),
    filter(unfiltered, !(PathwayType %in% subset.unfiltered$PathwayType)) %>% slice(1:topranked))

  # filter genes using study_universe
  genesets$universe <- genesets$Genes %>% unlist %>%
    intersect(., study_universe) %>% unique %>% length
  filtered <- genesets %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(Genes = list(intersect(Genes, study_universe)),
      n_geneset = length(Genes)) %>%
    dplyr::filter(n_geneset > min.set.size) %>%
    dplyr::mutate(n_input = length(genes),
      overlap = length(intersect(genes, Genes)),
      p_value = phyper(overlap - 1, n_geneset,
        universe - n_geneset, n_input, lower.tail = F)) %>%
    ungroup() %>%
    dplyr::group_by(PathwayType) %>%
    dplyr::mutate(bhfdr = p.adjust(p_value, method = "BH")) %>%
    arrange(PathwayType, p_value)
  subset.filtered <- summarise(filtered, sig=sum(bhfdr < alpha)) %>%
    filter(sig != 0) %>% dplyr::select(PathwayType)
  filtered <- rbind(filter(filtered, PathwayType %in% subset.filtered$PathwayType) %>% filter(bhfdr < alpha),
    filter(unfiltered, !(PathwayType %in% subset.filtered$PathwayType)) %>% slice(1:topranked))


  ## Proportion of mapped background set
  prop.mapped.backgroundSet <- unfiltered %>% summarise(int = round(100*length(intersect(unique(unlist(Genes)), study_universe))/length(study_universe), 1))

  prop <- cbind(mapped.genes.unfiltered, mapped.genes.filtered[,-1], prop.mapped.backgroundSet[,-1])
  colnames(prop) <- c("PathwayType", "mapped.genes.unfiltered", "mapped.genes.filtered", "prop.mapped.backgroundSet")

  return(list(unfiltered=unfiltered, filtered=filtered, prop=prop))
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
hypertest=function (genelist, genesets, cpus) {
  require(parallel)
  cl <- parallel::makeCluster(mc <- getOption("cl.cores", cpus))
  parallel::clusterExport(cl, varlist = c("genesets", "genelist"), envir = environment())
  pval <- parallel::parLapply(cl, genesets, function(genesetsi, genelist) {
    uni <- unlist(genesets)
    input <- genelist
    n_input = length(input)
    n_geneset = length(genesetsi)
    intersect = length(intersect(input, genesetsi))
    data.frame(pval=stats::phyper(intersect - 1, n_input, length(uni) - n_input,
      n_geneset, lower.tail = F),
      Genes=paste(intersect(input, genesetsi), collapse=";"),
      Overlap=paste0(intersect, "/", n_geneset))
  }, genelist)
  parallel::stopCluster(cl)
  pval %>% do.call(rbind, .) %>%
    as.data.frame %>% mutate(genesets = rownames(.),
      bhfdr = p.adjust(pval, "BH")) %>%
    arrange(pval) %>% dplyr::select(genesets,
      Overlap,
      pval, bhfdr, Genes)
}

