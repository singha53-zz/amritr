#' Convert .RCC file to data frame
#'
#' NanoString .RCC files to data frame
#' @param fileName file name including file path
#' @export
rccToDat = function(fileName) {
  library(dplyr); library(tidyr);
  lines <- data.frame(values = readLines(fileName))
  dat <- suppressWarnings(separate(data = lines, col = values,
    sep = ",", into = c("CodeClass", "Name", "Accession",
      "Count")))

  ind <- grep("<[A-Z]", dat$CodeClass)
  attr <- rep(NA, nrow(dat))
  for (i in 1:length(ind)) attr[ind[i]:nrow(dat)] <- grep("<[A-Z]",
    dat$CodeClass, value = TRUE)[i]
  dat <- dat %>% mutate(CodeClass = paste(CodeClass, gsub(" ",
    "", chartr("<>", "  ", attr)), sep = "_"), fileName = fileName)
  dat <- dat[-grep("<", dat$CodeClass), ]
  dat <- dat[!is.na(dat$Name), ]

  ## split flow cell data (properties) and biological (gene)
  ## data
  techDat <- dat[1:(grep("CodeClass", dat$CodeClass) - 1),
    ] %>% dplyr::select(-c(Accession:Count)) %>% spread(CodeClass,
      Name)
  bioDat <- dat[(grep("CodeClass", dat$CodeClass) + 1):nrow(dat),
    ]

  ## combine techDat and bioDat
  Dat <- full_join(techDat, bioDat, by = "fileName")

  return(Dat)
}
