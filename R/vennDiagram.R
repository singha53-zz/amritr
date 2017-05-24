#' Venn Diagram
#'
#' Venn diagram of 2, 3, or 4 sets
#' @param datList list of character vectors
#' @param circleNames names of list elements
#' @export
venndiagram = function(datList, circleNames){
  if(length(datList) != length(circleNames))
    stop("number of list elements does not match number of element names")

  if(length(datList) == 2){
    vennPlot <- vennDual(datList, circleNames)
  }
  if(length(datList) == 3){
    if(length(circleNames) != 3)
      stop("number of list elements does not match number of element names")
    vennPlot <- vennTriple(datList, circleNames)
  }
  if(length(datList) == 4){
    if(length(circleNames) != 4)
      stop("number of list elements does not match number of element names")
    vennPlot <- vennQuad(datList, circleNames)
  }

  ## determine specific overlaps
  combs <-
    unlist(lapply(1:length(datList),
      function(j) combn(names(datList), j, simplify = FALSE)),
      recursive = FALSE)
  names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
  elements <-
    lapply(combs, function(i) Setdiff(datList[i], datList[setdiff(names(datList), i)]))

  return(list(vennPlot=vennPlot, elements=elements))
}

#' Overlap between 2 sets
#'
#' @param datList list of character vectors
#' @export
vennDual = function(datList, circleNames){
  first <- datList[[1]]
  second <- datList[[2]]
  VennDiagram::draw.pairwise.venn(
    area1 = length(first),
    area2 = length(second),
    cross.area = length(intersect(first, second)),
    category = circleNames,
    fill = c("dodgerblue", "goldenrod1"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.pos = c(285, 105),
    cat.dist = 0.09,
    cat.just = list(c(-1, -1), c(1, 1)),
    ext.pos = 30,
    ext.dist = -0.05,
    ext.length = 0.85,
    ext.line.lwd = 2,
    ext.line.lty = "dashed"
  );
}

#' Overlap between 3 sets
#'
#' @param datList list of character vectors
#' @export
vennTriple = function(datList, circleNames){
  first <- datList[[1]]
  second <- datList[[2]]
  third <- datList[[3]]
  VennDiagram::draw.triple.venn(
    area1 = length(first),
    area2 = length(second),
    area3 = length(third),
    n12 = length(intersect(first, second)),
    n23 = length(intersect(second, third)),
    n13 = length(intersect(first, third)),
    n123 = length(Reduce(intersect, list(first, second, third))),
    category = circleNames,
    fill = c("dodgerblue", "goldenrod1", "darkorange1"),
    lty = "blank",
    cex = 2,
    cat.cex = 1,
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1")
  );
}

#' Overlap between 4 sets
#'
#' @param datList list of character vectors
#' @export
vennQuad = function(datList, circleNames){
  first <- datList[[1]]
  second <- datList[[2]]
  third <- datList[[3]]
  fourth <- datList[[4]]
  VennDiagram::draw.quad.venn(
    area1 = length(first),
    area2 = length(second),
    area3 = length(third),
    area4 = length(fourth),
    n12 = length(intersect(first, second)),
    n13 = length(intersect(first, third)),
    n14 = length(intersect(first, fourth)),
    n23 = length(intersect(second, third)),
    n24 = length(intersect(second, fourth)),
    n34 = length(intersect(third, fourth)),
    n123 = length(Reduce(intersect, list(first, second, third))),
    n124 = length(Reduce(intersect, list(first, second, fourth))),
    n134 = length(Reduce(intersect, list(first, third, fourth))),
    n234 = length(Reduce(intersect, list(second, third, fourth))),
    n1234 = length(Reduce(intersect, list(first, second, third, fourth))),
    category = circleNames,
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
    lty = "dashed",
    cex = 2,
    cat.cex = 0.75,
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3")
  );
}

#' Intersection function
#'
#' @param datList list of character vectors
#' @export
Intersect <- function (x) {
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

#' Union function
#'
#' @param datList list of character vectors
#' @export
Union <- function (x) {
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

#' Setdiff function
#'
#' @param datList list of character vectors
#' @export
Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's.
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}
