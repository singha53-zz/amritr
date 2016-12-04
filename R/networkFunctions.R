#' integrativePanels()
#'
#' build integrative panels (Concatenation, Ensemble, DIABLO)
#' @param X.train list of nxp datasets of length k
#' @param Y.train - n-vector of class labels
#' @param concat_alpha - 0 < value < 1
#' @param concat_lambda - value that controls the strength of the penalization
#' @param single_alphaList - list of alpha values of length k
#' @param single_lambdaList - list of lambda values of length k
#' @export
integrativePanels = function(X.train, Y.train, concat_alpha, concat_lambda, single_alphaList, single_lambdaList){
# Concatenation-Enet
combinedDat <- do.call(cbind, X.train)
result <- amritr::enet(X = combinedDat, Y = Y.train, alpha=concat_alpha, family="multinomial", lambda = concat_lambda, X.test = NULL,
  Y.test = NULL, filter = "none", topranked = 50)   # run once to determine lambda then set it
concat_enetPanel <- lapply(X.train, function(i){ intersect(colnames(i), result$enet.panel)})

# Ensemble-Enet
ensemble_enetPanel <- mapply(function(x, y, z){
  result <- amritr::enet(X = x, Y = Y.train, alpha=y, family="multinomial", lambda = z, X.test = NULL,
    Y.test = NULL, filter = "none", topranked = 50)
  result$enet.panel
}, x = X.train, y = single_alphaList, z = single_lambdaList)

# DIABLO-Enet
design <- matrix(1, nrow = length(X.train), ncol = length(X.train))
diag(design) <- 0
ncomp <- nlevels(Y.train) - 1
list.keepX = lapply(ensemble_enetPanel, function(i){
  rep(round(length(i)/ncomp, 0), ncomp)
})
diabloMod = mixOmics::block.splsda(X = X.train, Y = Y.train,
  ncomp = ncomp, keepX = list.keepX, design = design,
  scheme = "centroid")

diabloPanel <- lapply(diabloMod$loadings[-(length(X.train)+1)], function(x)
  unique(as.character(as.matrix(apply(x, 2, function(i) names(i)[which(i != 0)])))))

## Panels
panels = list(Concatenation = concat_enetPanel, Ensemble = ensemble_enetPanel,
  DIABLO = diabloPanel)

## Panel length
panelLength <- data.frame(Concatenation = unlist(lapply(concat_enetPanel, length)),
  Ensemble = unlist(lapply(ensemble_enetPanel, length)),
  DIABLO = unlist(lapply(diabloPanel, length)))

return(list(panels = panels, panelLength = panelLength, concat_lambda = result$lambda))
}

#' networkStats()
#'
#' compuate network statistics given a adjacency matrix
#' @param adjMat p x p matrix
#' @param mode = "graph"
#' @export
networkStats = function(adjMat, mode = "graph"){

nrelations <- network::network(adjMat,directed=FALSE)
#gplot(nrelations, usearrows=FALSE)

# Centrality
## Degree
nodeDegree <- sna::degree(nrelations)
# Betweenness
bet <- sna::betweenness(nrelations, gmode=mode) # Geographic betweenness
# Closeness
clo <- sna::closeness(nrelations) # Geographic closeness
# Closeness2
closeness2 <- function(x){ # Create an alternate closeness function!
  geo <- 1/sna::geodist(x)$gdist # Get the matrix of 1/geodesic distance
  diag(geo) <- 0 # Define self-ties as 0
  apply(geo, 1, sum) # Return sum(1/geodist) for each vertex
}
clo2 <- closeness2(nrelations) # Use our new function on contiguity data

# Geodesic distance
gdist <- sna::geodist(nrelations)$gdist # matrix of geodesic distances

# Eigenvector centrality score
eigenCentrality <- sna::evcent(nrelations)

# Harary graph centrality
# The Harary graph centrality of a vertex v is equal to 1/(max_u d(v,u)), where d(v,u)
# is the geodesic distance from v to u. Vertices with low graph centrality scores are
# likely to be near the “edge” of a graph, while those with high scores are likely to
# be near the “middle.” Compare this with closeness, which is based on the reciprocal
# of the sum of distances to all other vertices (rather than simply the maximum).
hararyCentrality <- sna::graphcent(nrelations)

# Prestige is the name collectively given to a range of centrality
# scores which focus on the extent to which one is nominated by
# others. default = domain: indegree within the reachability graph (Lin's unweighted measure)
prestigeScore <- sna::prestige(nrelations, cmode="domain")

# Stress Centrality
# The stress of a vertex, v, is given by
# C_S(v) = sum( g_ivj, i,j: i!=j,i!=v,j!=v)
# where g_ijk is the number of geodesics from i to k through j.
# Conceptually, high-stress vertices lie on a large number of
# shortest paths between other vertices; they can thus be thought
# of as “bridges” or “boundary spanners.” Compare this with
# betweenness, which weights shortest paths by the inverse of
# their redundancy.
stressScore <- sna::stresscent(nrelations)     #Compute stress scores

# Graph-level indicies
graphDensity <- sna::gden(nrelations) # graph density
## Reciprocity
reciprocityMeasure <- c("dyadic", "dyadic.nonnull", "edgewise", "edgewise.lrr", "correlation")
recipScore <- unlist(lapply(reciprocityMeasure, function(i){
  sna::grecip(nrelations, measure = i)
}))
names(recipScore) <- reciprocityMeasure

transitivityMeasure <- c("weak", "strong", "weakcensus", "strongcensus", "rank", "correlation")
transScore <- unlist(lapply(transitivityMeasure, function(i){
  sna::gtrans(nrelations, mode = mode, measure = i)
}))
names(transScore) <- transitivityMeasure

# Census
dyadCensus <- sna::dyad.census(nrelations) # M,A,N counts
triadCensus <- sna::triad.census(nrelations, mode = mode) # Directed triad census

# Compute k-path or k-cycle census statistics
#pathCensus <- kpath.census(nrelations, mode = mode, maxlen=maxlen, tabulate.by.vertex=FALSE) # Count paths of length <=6
#cycleCensus <- kcycle.census(nrelations, mode = mode, maxlen=maxlen, tabulate.by.vertex=FALSE) # Count cycles of length <=6

# Cycle Census information
# A (maximal) clique is a maximal set of mutually adjacenct
# vertices. Cliques are important for their role as cohesive
# subgroups
cliques <- sna::clique.census(nrelations, mode = mode, tabulate.by.vertex=FALSE, enumerate=FALSE)$clique.count # Find maximal cliques

# isolate: An isolated vertex is a vertex with
# degree zero; that is, a vertex that is not an
# endpoint of any edge
isol <- sna::isolates(nrelations) # Get the entire list of isolates
isolatedVertices <- network::network.vertex.names(nrelations)[isol]

# Compure k-Core
kc <- sna::kcores(nrelations, mode = mode, cmode="indegree")

nodeIndices <- rbind(nodeDegree, bet, clo, clo2, eigenCentrality, hararyCentrality, prestigeScore, stressScore, kc)

return(list(nodeIndices = nodeIndices, gdist = gdist, graphDensity = graphDensity, recipScore = recipScore,
  transScore = transScore, dyadCensus = dyadCensus, triadCensus = triadCensus, cliques = cliques, isolatedVertices = isolatedVertices))
}

#' graphIndices()
#'
#' build integrative panels (Concatenation, Ensemble, DIABLO)
#' @param panels - list of panels (character vector corresponding to column names of X.train)
#' @param X.train - list of datasets (n x p)
#' @param cut-off - of person correlatio
#' @param concat_lambda - value that controls the strength of the penalization
#' @param single_alphaList - list of alpha values of length k
#' @param single_lambdaList - list of lambda values of length k
#' @export
graphIndices = function(panels = panels, X.train = X.train, cutoff = cutoff){
  library(dplyr)
  ## Determine adjacency matrices
  adjMat <- lapply(panels, function(i){
    adjMat = cor(do.call(cbind, mapply(function(x, y){
      y[, x]
    }, x = i, y = X.train)))
    adjMat[abs(adjMat) < cutoff] <- 0
    adjMat[abs(adjMat) > cutoff] <- 1
    adjMat
  })

  # Estimate graph statistics
  graphs <- lapply(adjMat, function(x) amritr::networkStats(adjMat = x))

  ## graphDensity, recipScore, transScore, cliques
  gIndices <- lapply(graphs, function(i){
    data.frame(Value = c(i$graphDensity, i$recipScore["edgewise.lrr"], i$transScore["weakcensus"], length(i$cliques)),
      Statistic = c("Graph Density", "Edgewise.lrr", "Transitivity", "NumOfCliques"))

  })
  gIndicesDat <- do.call(rbind, gIndices) %>% mutate(Method = rep(names(gIndices), each = 4))

  # Dyads and Triads
  dyads <- do.call(rbind, lapply(graphs, function(x){x$dyadCensus}))
  dyads <- as.data.frame(dyads)[,-2] %>% mutate(Method = names(graphs)) %>%
    tidyr::gather(Type, Number, -Method)
  triads <- do.call(rbind, lapply(graphs, function(x){x$triadCensus}))
  triads <- as.data.frame(triads) %>% mutate(Method = names(graphs)) %>%
    tidyr::gather(Type, Number, -Method)

  # Number of isolates
  isolatedFeat <- lapply(graphs, function(x){ x$isolatedVertices})
  isolates <- do.call(rbind, lapply(isolatedFeat, function(x){
    unlist(lapply(X.train, function(y){
      length(intersect(x, colnames(y)))
    }))
  })) %>% as.data.frame %>% mutate(Method = names(graphs)) %>%
    tidyr::gather(Dataset, NumOfIsolates, -Method)

  return(list(adjMat = adjMat, gIndicesDat = gIndicesDat, dyads = dyads, triads = triads, isolates = isolates))
}
