## functions implementing the GMA-SMA models from Forst et al (2019)

## ============================================================================
## code to handle input data sets

##' Reformat composition of mixtures into a list
##'
##' Reformats composition of mixtures into a list.
##' @param mixtures vector, matrix, data.frame or list containing of mixtures
##' @param sep required only if \code{mixtures} is a vector; will be given to \code{strsplit}
##' @return list of vectors, with one component per mixture, the vector containing the mixture components
##' @author Timothee Flutre
##' @seealso \code{\link{getMixturesPerGeno}}
##' @examples
##' nbGenos <- 25
##' levGenos <- sprintf(fmt=paste0("geno%0", floor(log10(nbGenos))+1, "i"),
##'                     1:nbGenos)
##' nbMixes <- 75 # only binary and balanced
##' design <- getDesignBinaryVarMix(levGenos, nbMixes, seed=12345)
##' mixtures <- getMixtureList(design$combs)
##' str(mixtures, list.len=10)
##' @export
getMixtureList <- function(mixtures, sep = NULL) {
  out <- list()

  if (is.factor(mixtures))
    mixtures <- as.character(mixtures)

  if (is.vector(mixtures) & !is.list(mixtures)) {
    stopifnot(!is.null(sep))
    mixtures <- as.character(mixtures)
    out <- strsplit(mixtures, paste0("\\", sep))
    if (!is.null(names(mixtures))) {
      names(out) <- names(mixtures)
    } else {
      names(out) <- mixtures
    }
  } else {
    if (is.matrix(mixtures)) {
      stopifnot(!is.null(rownames(mixtures)))
      mixtures <- as.data.frame(mixtures)
    }
    if (is.data.frame(mixtures)) {
      colnames(mixtures) <- NULL
      stopifnot(!is.null(rownames(mixtures)))
      mix_names <- rownames(mixtures)
      out <- lapply(1:nrow(mixtures), function(i) {
        unlist(mixtures[i, ])
      })
      names(out) <- mix_names
    } else if (is.list(mixtures)) {
      if (all(sapply(mixtures, is.factor))) {
        mixtures <- lapply(mixtures, as.character)
      }
      stopifnot(
        all(sapply(mixtures, is.vector)),
        all(!sapply(mixtures, anyDuplicated))
      )
      out <- mixtures
    }
  }

  out <- lapply(out, as.character)

  return(out)
}

##' Get mixture(s) per genotype
##'
##' Returns the list of mixtures per genotype.
##' @param mix2genos list with one component per mixture, listing the genotypes it contains (output of \code{\link{getMixtureList}})
##' @return list with one vector per genotype
##' @author Timothee Flutre
##' @seealso \code{\link{getMixtureList}}
##' @examples
##' nbGenos <- 25
##' levGenos <- sprintf(fmt=paste0("geno%0", floor(log10(nbGenos))+1, "i"),
##'                     1:nbGenos)
##' nbMixes <- 75 # only binary and balanced
##' design <- getDesignBinaryVarMix(levGenos, nbMixes, seed=12345)
##' mixtures <- getMixtureList(design$combs)
##' str(mixtures, list.len=10)
##' geno2mixes <- getMixturesPerGeno(mixtures)
##' geno2mixes[c(1,2)]
##' table(sapply(geno2mixes, length))
##' @export
getMixturesPerGeno <- function(mix2genos) {
  stopifnot(is.list(mix2genos))
  genos <- unique(sort(do.call(c, mix2genos)))
  mix2genos <- lapply(mix2genos, paste, collapse = ",")
  vec_mix2genos <- do.call(c, mix2genos)
  out <- Map(function(geno) {
    idx <- grep(geno, vec_mix2genos)
    names(mix2genos)[idx]
  }, genos)
  return(out)
}

##' Pivot mixture data into the long format
##'
##' Pivots mixture data into the long format, from one row per stand into one row per component of each stand (i.e., several rows if the stand is a mixture).
##' @param df data frame
##' @param genos named list of vectors, one per species, containing all the identifiers of the genotypes of the given species
##' @param colC name of the column containing the component(s) of each stand
##' @param sep separator used to separate the genotype names in \code{"col"}; see the \code{"split"} argument of \code{\link{strsplit}}
##' @param prefixY prefix of the columns containing the yield data (one column per mixture component); it is assumed that the components in \code{colC} are in the same order as the yield data
##' @param sepY separator used in the name of the new "yield" columns made by concatenating \code{colY} and \code{colIDfocal} or \code{colIDneighbors}
##' @param colIDfocal name of the column in the output that will contain the identifiers of the focal genotypes
##' @param colIDneighbors name of the column in the output that will contain the identifiers of the neighbor genotypes; for monovarietal stands, the entries in this neighbor column will be identical to the entries in the focal column
##' @return data.frame with more rows
##' @seealso \code{\link{pivotMixData2Wide}}
##' @examples
##' ## example of species mixtures:
##' (dat <- data.frame(name=c("wheat01-pea02","wheat01-pea01", "wheat01"),
##'                    stand=c("inter","inter","sole"),
##'                    x=c(1, 1, 1),
##'                    y=c(1, 2, 3),
##'                    "yield_cereal" = c(10, 11, 20),
##'                    "yield_legume" = c(9, 8, NA),
##'                    check.names = FALSE,
##'                    stringsAsFactors=TRUE))
##' pivotMixData2Long(df=dat, colC="name", sep="-", prefixY="yield",
##'                   genos=list("cereal"=c("wheat01"),
##'                              "legume"=c("pea01","pea02")))
##'
##' ## example of varietal mixtures:
##' (dat <- data.frame(name=c("g1+g2","g1+g2","g1","g2"),
##'                    "yield_focal" = c(10, 12, 20, 18),
##'                    "yield_neighbor" = c(11, 9, NA, NA),
##'                    check.names = FALSE,
##'                    stringsAsFactors=TRUE))
##' pivotMixData2Long(df=dat, colC="name", sep="+", prefixY="yield",
##'                   genos=list(c("g1","g2")))
##' @export
pivotMixData2Long <- function(df, genos, colC, sep = "-",
                              prefixY = NULL, sepY = "_",
                              colIDfocal = "focal", colIDneighbors = "neighbor") {
  if (FALSE) { # debug
    df <- dat
    genos <- list("cereal" = c("wheat01"), "legume" = c("pea01", "pea02"))
    colC <- "name"
    sep <- "-"
    prefixY <- "yield"
    colIDfocal <- "focal"
    colIDneighbors <- "neighbor"
  }
  stopifnot(
    is.data.frame(df),
    colC %in% colnames(df),
    all(!c(colIDfocal, colIDneighbors) %in% colnames(df))
  )
  colY <- NA
  if (!is.null(prefixY)) {
    if (length(genos) == 1) {
      colY <- grep(paste0("^", prefixY), colnames(df), value = TRUE)
    } else {
      colY <- paste0(prefixY, sepY, names(genos))
    }
    stopifnot(colY %in% colnames(df))
  }

  comps <- strsplit(as.character(df[[colC]]), split = paste0("\\", sep))
  if (!is.null(prefixY)) {
    yields <- apply(df[, colY], 1, as.vector, simplify = FALSE)
    yields <- lapply(yields, function(x) {
      x[which(!is.na(x))]
    })
    stopifnot(all(sapply(comps, length) == sapply(yields, length)))
    yields <- Map(function(x1, x2) {
      setNames(x1, x2)
    }, yields, comps)
  }

  out <- lapply(seq_along(comps), function(i) {
    n <- length(comps[[i]])
    tmp <- as.data.frame(matrix(nrow = n, ncol = 2))
    colnames(tmp) <- c(colIDfocal, colIDneighbors)
    tmp[[colIDfocal]] <- comps[[i]]
    for (j in 1:n) {
      tmp[[colIDneighbors]][j] <- paste0(comps[[i]][comps[[i]] != comps[[i]][j]], collapse = sep)
    }
    cbind(
      tmp,
      df[rep(i, n), ]
    )
  })

  out <- do.call(rbind, out)

  isMono <- (out[[colIDneighbors]] == "")
  if (any(isMono)) {
    idx <- which(isMono)
    out[[colIDneighbors]][idx] <- out[[colIDfocal]][idx]
  }

  if (!is.null(prefixY)) {
    out[[prefixY]] <- do.call(c, yields)
    out[colY] <- NULL
  }

  out[[colIDfocal]] <- factor(out[[colIDfocal]])
  out[[colIDneighbors]] <- factor(out[[colIDneighbors]])
  rownames(out) <- NULL

  return(out)
}

##' Pivot mixture data into the wide format
##'
##' Pivots mixture data into the wide format, from one row per component of each stand (i.e., several rows if the stand is a mixture) into one row per stand.
##' @param df data frame
##' @param colIDstand name of the column containing the identifier of each stand
##' @param colIDfocal name of the column containing the identifiers of the focal genotypes
##' @param colIDneighbors optional name of the column containing the identifiers of the neighbor genotypes
##' @param colPlot name of the column(s) allowing to uniquely identify a plot
##' @param colY name of the column containing the yield data
##' @param sepY separator used in the name of the new "yield" columns made by concatenating \code{colY} and \code{colIDfocal} or \code{colIDneighbors}
##' @param sepFocalNeighbors in case \code{colIDneighbors} is unspecified or does not exist in \code{df}, the distinction between focal and neighbor(s) will be retrieved by splitting \code{colIDstand} using \code{sepFocalNeighbors}
##' @return data.frame
##' @seealso \code{\link{pivotMixData2Long}}
##' @author Timothee Flutre
##' @examples
##' ## only binary mixtures:
##' dat0 <- data.frame(
##'   focal = c("wheat01", "pea02", "wheat01", "pea01"),
##'   neighbor = c("pea02", "wheat01", "pea01", "wheat01"),
##'   name = c(
##'     rep("wheat01-pea02", 2),
##'     rep("wheat01-pea01", 2)
##'   ),
##'   stand = rep("inter", 4),
##'   x = 1,
##'   y = c(1, 1, 2, 2),
##'   yield = c(10, 11, 12, 13),
##'   stringsAsFactors = TRUE
##' )
##' dat0
##' (dat1 <- pivotMixData2Wide(dat0, colIDstand="name"))
##'
##' ## binary mixtures and a monovarietal:
##' dat0 <- data.frame(
##'   focal = c("wheat01", "pea02", "wheat01", "pea01", "wheat01"),
##'   neighbor = c("pea02", "wheat01", "pea01", "wheat01", "wheat01"),
##'   name = c(
##'     rep("wheat01-pea02", 2),
##'     rep("wheat01-pea01", 2),
##'     "wheat01"
##'   ),
##'   stand = c(rep("inter", 4), "sole"),
##'   x = 1,
##'   y = c(1, 1, 2, 2, 3),
##'   yield = c(10, 11, 12, 13, 20.5),
##'   stringsAsFactors = TRUE
##' )
##' dat0
##' (dat1 <- pivotMixData2Wide(dat0, colIDstand="name"))
##'
##' ## reverse conversion
##' dat1v2 <- dat1
##' colnames(dat1v2)[colnames(dat1v2) == "yield_focal"] <- "yield-cereal"
##' colnames(dat1v2)[colnames(dat1v2) == "yield_neighbor"] <- "yield-legume"
##' dat1v2$yield <- apply(dat1v2[,c(7,8)], 1, sum, na.rm=TRUE)
##' (dat1v2 <- dat1v2[,-c(1,2)])
##' (dat2 <- pivotMixData2Long(dat1v2, colC="name", sepY="-",
##'                            genos=list("cereal"=c("wheat01","wheat02"),
##'                                       "legume"=c("pea01","pea02"))))
##' all.equal(dat2, dat0)
##' @export
pivotMixData2Wide <- function(df, colIDstand = "ID",
                              colIDfocal = "focal",
                              colIDneighbors = "neighbor",
                              colPlot = c("x", "y"),
                              colY = "yield", sepY = "_",
                              sepFocalNeighbors = NULL) {
  if (FALSE) { # debug
    df <- dat
    colIDstand <- "name"
    colIDfocal <- "focal"
    colIDneighbors <- "neighbor"
    colPlot <- c("x", "y")
    colY <- "yield"
  }
  stopifnot(
    is.data.frame(df),
    all(c(colIDstand, colIDfocal, colPlot, colY) %in% colnames(df)),
    all(!is.na(df[[colIDstand]])),
    all(!is.na(df[[colIDfocal]])),
    all(!is.na(df[[colIDneighbors]]))
  )

  ## add "neighbors" column if needed:
  if (any(is.null(colIDneighbors),
          ! colIDneighbors %in% colnames(df))) {
    stopifnot(! is.null(sepFocalNeighbors))
    if (is.null(colIDneighbors))
      colIDneighbors <- "neighbor"
    stopifnot(! colIDneighbors %in% colnames(df))
    mix2genos <- getMixtureList(df[[colIDstand]], sep = sepFocalNeighbors)
    neighbors <- Map(function(i){
      mix2genos[[i]][which(mix2genos[[i]] != df[[colIDfocal]][i])]
    }, 1:nrow(df))
    df[[colIDneighbors]] <- sapply(neighbors, paste, collapse = sepFocalNeighbors)
    df[[colIDneighbors]] <- factor(df[[colIDneighbors]])
  }

  ## define keys (uniquely identifying plots):
  stopifnot(!"key" %in% colnames(df))
  df[["key"]] <- apply(df[, colPlot, drop = FALSE], 1, paste0, collapse = "_")
  df[["key"]] <- factor(df[["key"]])

  ## find max mixture order (nb of comps per mix):
  tmp <- Map(function(k) {
    df[which(df$key == k), ]
  }, levels(df[["key"]]))
  max_mix_order <- max(sapply(tmp, nrow))

  ## set yield column names:
  colYnew <- paste0(colY, sepY, colIDfocal)
  if (max_mix_order == 2) {
    colYnew <- c(
      colYnew,
      paste0(colY, sepY, colIDneighbors)
    )
  } else {
    colYnew <- c(
      colYnew,
      paste0(colY, sepY, colIDneighbors, sepY, 1:(max_mix_order - 1))
    )
  }
  df[, colYnew] <- NA

  ## make the output:
  out <- Map(function(k) {
    tmp <- df[which(df$key == k), ]
    idx <- 1:nrow(tmp) # to handle stands of any order (mono, binary mix, ternary mix, etc)
    out_k <- tmp[1, colnames(df)[colnames(df) != colY]]
    out_k[, colYnew[idx]] <- tmp[[colY]][idx]
    out_k
  }, levels(df[["key"]]))
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[["key"]] <- NULL
  out <- droplevels(out)

  return(out)
}

## ============================================================================
## code to make design matrices

##' @noRd
getNeighbors <- function(df, col, sep) {
  neighbors <- strsplit(as.character(df[[col]]), paste0("\\", sep))
  names(neighbors) <- as.character(df[[col]])
  neighbors <- lapply(neighbors, unique)

  if (all(neighbors[[1]] == "")) {
    msg <- paste0(
      "split of column '", col, "' did not work,",
      " check that your sep ('", sep, "') works with strsplit"
    )
    stop(msg)
  }

  return(neighbors)
}

##' @noRd
getGenos <- function(neighbors) {
  return(sort(unique(do.call(c, neighbors))))
}

##' Design matrix of GMAs for GMA-SMA models
##'
##' Makes a design matrix with the contribution of general mixing abilities (GMAs).
##' The design is based on the names of pure and mixed stands, in which the genotype names are separated by a specific symbol.
##' @param df data frame
##' @param col name of the column containing the names of pure and mixed stands
##' @param sep separator used to separate the genotype names in \code{"col"}
##' @return matrix
##' @seealso \code{\link{mkZGMA}}, \code{\link{fitGMASMA}}
##' @author Timothee Flutre
##' @examples
##' dat <- data.frame(mix=paste0("mix", 1:3),
##'                   varieties=c("var3","var1-var3","var1-var2"),
##'                   pheno=c(10, 7, 9))
##' (Z_GMA <- mkZGMA(df=dat, col="varieties", sep="-"))
##' @export
mkZGMA <- function(df, col, sep = ",") {
  stopifnot(
    is.data.frame(df),
    col %in% colnames(df)
  )

  neighbors <- getNeighbors(df, col, sep)
  genos <- getGenos(neighbors)

  Z_GMA <- matrix(0, nrow = nrow(df), ncol = length(genos))
  colnames(Z_GMA) <- genos

  idx <- which(sapply(neighbors, length) == 1) # pure stands
  for (i in idx) {
    Z_GMA[i, neighbors[[i]][1]] <- 1
  }

  idx <- which(sapply(neighbors, length) > 1) # mixed stands
  for (i in idx) {
    Z_GMA[i, neighbors[[i]]] <- 1 / length(neighbors[[i]])
  }

  stopifnot(all(rowSums(Z_GMA) == 1))

  return(Z_GMA)
}

##' Design matrix of SMAs for GMA-SMA models
##'
##' Makes a design matrix with the contribution of specific mixing abilities (SMAs).
##' Before Forst et al (2019), there was only one kind of SMA, the inter-genotypic SMA (SMA_ij).
##' Forst et al (2019) introduced a second kind of SMA, the intra-genotypic SMA (SMA_ii).
##' As a result, various design matrices can be made depending on which kinds of SMAs are modeled.
##' The design is based on the names of pure and mixed stands, in which the genotype names are separated by a specific symbol.
##' @param df data frame
##' @param col name of the column containing the names of pure and mixed stands
##' @param sep separator used to separate the genotype names in \code{"col"}
##' @param inc_SMA_ii specify how the intra-genotypic SMA should be included, with "no" meaning no SMA_ii, "only_pur" meaning that the SMA_ii is only included in the pure stands (this is model 2 of Forst et al, 2019), "pur_mix" meaning that the SMA_ii is included in both the pure and the mixed stands (this is model 3 of Forst et al, 2019)
##' @param skipUnusedCols skip unused columns, i.e., full of zeroes and not present in the data frame
##' @return matrix
##' @seealso \code{\link{mkZGMA}}, \code{\link{fitGMASMA}}
##' @author Timothee Flutre
##' @examples
##' dat <- data.frame(mix=paste0("mix", 1:3),
##'                   varieties=c("var3","var1-var3","var1-var2"),
##'                   pheno=c(10, 7, 9))
##' (Z_SMA <- mkZSMA(df=dat, col="varieties", sep="-", inc_SMA_ii="only_pur")) # model 2
##' (Z_SMA <- mkZSMA(df=dat, col="varieties", sep="-", inc_SMA_ii="pur_mix")) # model 3
##' @export
mkZSMA <- function(df, col, sep = ",", inc_SMA_ii = "no", skipUnusedCols = TRUE) {
  stopifnot(
    is.data.frame(df),
    col %in% colnames(df),
    inc_SMA_ii %in% c("no", "only_pur", "pur_mix"),
    is.logical(skipUnusedCols)
  )

  neighbors <- getNeighbors(df, col, sep)
  genos <- getGenos(neighbors)

  colNames <- expand.grid(genos, genos)
  colNames <- apply(colNames, 1, function(x) {
    paste(sort(x), collapse = sep)
  })
  colNames <- colNames[!duplicated(colNames)]

  Z_SMA <- matrix(0, nrow = nrow(df), ncol = length(colNames))
  colnames(Z_SMA) <- colNames

  ## "classical" model -> model 2'
  if (inc_SMA_ii == "no") {
    allPairs <- lapply(neighbors, function(x) {
      if (length(x) == 1) {
        paste(rep(x, 2), collapse = sep)
      } else {
        tmp <- t(utils::combn(sort(x), 2))
        paste(tmp[, 1], tmp[, 2], sep = sep)
      }
    })
    for (i in 1:nrow(Z_SMA)) {
      Z_SMA[i, allPairs[[i]]] <- 1 / length(allPairs[[i]])
    }
    isPure <- (sapply(
      strsplit(colnames(Z_SMA), paste0("\\", sep)),
      anyDuplicated
    ) == 2)
    colsToKeep <- (!isPure) # no SMA_ii in model 2
    Z_SMA <- Z_SMA[, colsToKeep]
  }

  ## model 2 of Forst et al (2019)
  if (inc_SMA_ii == "only_pur") {
    allPairs <- lapply(neighbors, function(x) {
      if (length(x) == 1) {
        paste(rep(x, 2), collapse = sep)
      } else {
        tmp <- t(utils::combn(sort(x), 2))
        paste(tmp[, 1], tmp[, 2], sep = sep)
      }
    })
    for (i in 1:nrow(Z_SMA)) {
      Z_SMA[i, allPairs[[i]]] <- 1 / length(allPairs[[i]])
    }
    isPure <- (sapply(
      strsplit(colnames(Z_SMA), paste0("\\", sep)),
      anyDuplicated
    ) == 2)
    isEmpty <- (colSums(Z_SMA) == 0)
    colsToRmv <- (isPure & isEmpty) # rm SMA_ii in model 2 if empty
    if (any(colsToRmv)) {
      Z_SMA <- Z_SMA[, -which(colsToRmv)]
    }
  }

  ## model 3 of Forst et al (2019)
  if (inc_SMA_ii == "pur_mix") {
    allPairs <- lapply(neighbors, function(x) {
      tmp <- expand.grid(x, x)
      isPure <- (apply(tmp, 1, anyDuplicated) != 0)
      tmp <- apply(tmp, 1, function(x) {
        paste(sort(x), collapse = sep)
      })
      isUniq <- (!duplicated(tmp))
      data.frame(
        id = tmp[isUniq],
        isPure = isPure[isUniq]
      )
    })
    for (i in 1:nrow(Z_SMA)) {
      if (all(allPairs[[i]]$isPure)) {
        Z_SMA[i, allPairs[[i]]$id[allPairs[[i]]$isPure]] <- 1
      } else {
        Z_SMA[i, allPairs[[i]]$id[allPairs[[i]]$isPure]] <- 1 / length(neighbors[[i]])^2
        Z_SMA[i, allPairs[[i]]$id[!allPairs[[i]]$isPure]] <- 2 / length(neighbors[[i]])^2
      }
    }
  }

  if (skipUnusedCols) {
    isEmpty <- (colSums(Z_SMA) == 0)
    isPure <- (sapply(
      strsplit(colnames(Z_SMA), paste0("\\", sep)),
      anyDuplicated
    ) == 2)
    colsToKeep <- (!isEmpty) | (isPure) # keep SMA_ii in model 3
    Z_SMA <- Z_SMA[, colsToKeep]
  }

  stopifnot(all(rowSums(Z_SMA) %in% c(0, 1)))

  return(Z_SMA)
}

##' Design matrices of SMAs for GMA-SMA models
##'
##' Makes the design matrices with the contribution of specific mixing abilities (SMAs).
##' More details in the help of \code{\link{mkZSMA}} and in the vignette.
##' @param df data frame
##' @param col name of the column containing the names of pure and mixed stands
##' @param sep separator used to separate the genotype names in \code{"col"}
##' @param models vector of model identifiers; see the vignette
##' @param verbose verbosity level
##' @return list of SMA design matrices
##' @author Timothee Flutre
##' @examples
##' dat <- data.frame(mix=paste0("mix", 1:3),
##'                   varieties=c("var3","var1-var3","var1-var2"),
##'                   pheno=c(10, 7, 9))
##' (listZSMA <- mkAllZSMA(df=dat, col="varieties", sep="-"))
##' @export
mkAllZSMA <- function(df, col, sep,
                      models = c("2p", "2", "2pp", "3", "3p"),
                      verbose = FALSE) {
  stopifnot(all(models %in% c("2p", "2", "2pp", "3", "3p")))

  out <- list()

  for (model in models) {
    if (model == "2p") {
      if (verbose) {
        print("model 2'")
      }
      out$SMA_mod2p <- mkZSMA(df, col, sep, inc_SMA_ii = "no")
    } else if (model %in% c("2", "2pp")) {
      if (!"SMA_mod2" %in% names(out)) {
        if (verbose) {
          print("model 2")
        }
        out$SMA_mod2 <- mkZSMA(df, col, sep, inc_SMA_ii = "only_pur")
      }

      if (model == "2pp") {
        if (verbose) {
          print("model 2''")
        }
        isPure <- (sapply(
          strsplit(colnames(out$SMA_mod2), paste0("\\", sep)),
          anyDuplicated
        ) == 2)
        out$SMA_mod2pp_ij <- out$SMA_mod2[, !isPure]
        out$SMA_mod2pp_ii <- out$SMA_mod2[, isPure]
      }
    } else if (model %in% c("3", "3p")) {
      if (!"SMA_mod3" %in% names(out)) {
        if (verbose) {
          print("model 3")
        }
        out$SMA_mod3 <- mkZSMA(df, col, sep, inc_SMA_ii = "pur_mix")
      }

      if (model == "3p") {
        if (verbose) {
          print("model 3'")
        }
        isPure <- (sapply(
          strsplit(colnames(out$SMA_mod3), paste0("\\", sep)),
          anyDuplicated
        ) == 2)
        out$SMA_mod3p_ij <- out$SMA_mod3[, !isPure]
        out$SMA_mod3p_ii <- out$SMA_mod3[, isPure]
      }
    }
  }

  return(out)
}

##' Design matrices for DBV-SBV models in an inter-specific trial
##'
##' Makes the design matrices for direct and social breeding values per species when the trial includes binary intercrops and, possibly, sole crops.
##' @param df data frame with one row per sole-crop plot or two rows per binary intercrop plot, with a factor column corresponding to the focal identifiers ("DBV" a.k.a. "DGE" or "Pr") and a column corresponding to the "neighbor(s)" identifiers ("SBV" a.k.a. "IGE" or "As")
##' @param genosPerSp list of two components, one per species, each being a vector with the genotype identifiers of the given species (they will be sorted)
##' @param colIDfocal name of the column containing the identifiers of the focal genotypes
##' @param colIDneighbors name of the column containing the identifiers of the neighbour genotypes
##' @return list of two components, one per species, each being a list of two matrices, the first for the DBVs and the second for the SBVs; the column names of these matrices correspond to the genotype identifiers after sorting
##' @author Timothee Flutre
##' @examples
##' (dat <- data.frame(focal=c("wheat01", "pea02",
##'                            "wheat01", "pea01",
##'                            "wheat01"),
##'                    neighbors=c("pea02", "wheat01",
##'                                "pea01", "wheat01",
##'                                "wheat01"),
##'                    name=c("wheat01-pea02", "wheat01-pea02",
##'                           "wheat01-pea01", "wheat01-pea01",
##'                           "wheat01-wheat01"),
##'                    stand = c("inter", "inter", "inter", "inter", "sole"),
##'                    species = c("wheat", "pea", "wheat", "pea", "wheat")))
##' (out <- mkZinterspe(dat,
##'                     list("wheat"="wheat01",
##'                          "pea"=c("pea01","pea02"))))
##' @export
mkZinterspe <- function(df, genosPerSp, colIDfocal = "focal", colIDneighbors = "neighbors") {
  stopifnot(
    is.data.frame(df),
    is.list(genosPerSp),
    length(genosPerSp) >= 2,
    all(c(colIDfocal, colIDneighbors) %in% colnames(df))
  )

  out <- list()

  genosPerSp <- lapply(genosPerSp, sort)

  ## check if there are missing genotypes per species
  ## by partitioning the row indices into focal/neighbor,
  ## and identify sole/inter-crops per species, too
  focsIdxPerSp <- list()
  neisIdxPerSp <- list()
  solesIdxPerSp <- list()
  intersIdxPerSp <- list()
  for (sp in names(genosPerSp)) {
    idxF <- which(df[[colIDfocal]] %in% genosPerSp[[sp]])
    focsIdxPerSp[[sp]] <- idxF
    genosF <- unique(sort(as.character(df[[colIDfocal]][idxF])))
    missGenosF <- genosF[!genosF %in% genosPerSp[[sp]]]
    if (length(missGenosF) > 0) {
      msg <- paste0(
        length(missGenosF), " focal genotype",
        ifelse(length(missGenosF) > 1, "s", ""),
        " from species ", sp, " missing from",
        " argument genosPerSp"
      )
      stop(msg)
    }
    idxN <- which(df[[colIDneighbors]] %in% genosPerSp[[sp]])
    neisIdxPerSp[[sp]] <- idxN
    genosN <- unique(sort(as.character(df[[colIDneighbors]][idxN])))
    missGenosN <- genosN[!genosN %in% genosPerSp[[sp]]]
    if (length(missGenosN) > 0) {
      msg <- paste0(
        length(missGenosN), " neighbor genotype",
        ifelse(length(missGenosN) > 1, "s", ""),
        " from species ", sp, " missing from",
        " argument genosPerSp"
      )
      stop(msg)
    }
    idxS <- which(df[[colIDfocal]] %in% genosPerSp[[sp]] &
      df[[colIDneighbors]] == df[[colIDfocal]])
    solesIdxPerSp[[sp]] <- idxS
    idxI <- which(df[[colIDfocal]] %in% genosPerSp[[sp]] &
      df[[colIDneighbors]] != df[[colIDfocal]])
    intersIdxPerSp[[sp]] <- idxI
  }

  ## make the design matrices for all species jointly
  Z_F_all <- model.matrix(formula(paste0("~ 0 + ", colIDfocal)), data = df)
  colnames(Z_F_all) <- lapply(colnames(Z_F_all), function(x) gsub(colIDfocal, "", x))
  Z_N_all <- model.matrix(formula(paste0("~ 0 + ", colIDneighbors)), data = df)
  colnames(Z_N_all) <- lapply(colnames(Z_N_all), function(x) gsub(colIDneighbors, "", x))

  ## make the design matrices per species
  for (sp in names(genosPerSp)) {
    out[[sp]] <- list()
    Z_F <- Z_F_all[, genosPerSp[[sp]], drop = FALSE]
    out[[sp]]$Z_DBV <- Z_F
    Z_N <- Z_N_all[, genosPerSp[[sp]], drop = FALSE]
    out[[sp]]$Z_SBV <- Z_N
  }

  ## set neighbors at 0 for sole crops
  for (sp in names(solesIdxPerSp)) {
    idxS <- solesIdxPerSp[[sp]]
    if (length(idxS) == 0) {
      next
    }
    out[[sp]]$Z_SBV[idxS, ] <- 0
  }

  return(out)
}

##' Design matrix for the SIS component in an inter-specific trial
##'
##' Makes the "SIS" design matrix (social interspecific effects) when the trial includes binary intercrops and sole crops.
##' @param df data frame with one row per sole-crop plot or two rows per binary intercrop plot, with a factor column corresponding to the focal identifiers ("DBV" a.k.a. "DGE" or "Pr")
##' @param colS name of the column containing the species of the focal genotype
##' @param species the species for which the SIS design matrix should be made
##' @param weight weight to be used for the mixed stands (by default, "1" is used for the pure stands)
##' @param colIDfocal name of the column containing the names for the focal genotype
##' @param colStand name of the column indicating the stand
##' @param levOneComp name indicating that the stand has a single component (usually "pure" or "sole")
##' @param levSevComps name indicating that the stand has several components (usually "mixed" or "inter")
##' @return matrix
##' @author Timothee Flutre
##' @examples
##' (dat <- data.frame(focal=c("wheat01", "pea02",
##'                            "wheat01", "pea01",
##'                            "wheat01"),
##'                    neighbors=c("pea02", "wheat01",
##'                                "pea01", "wheat01",
##'                                "wheat01"),
##'                    name=c("wheat01-pea02", "wheat01-pea02",
##'                           "wheat01-pea01", "wheat01-pea01",
##'                           "wheat01-wheat01"),
##'                    stand = c("inter", "inter", "inter", "inter", "sole"),
##'                    species = c("wheat", "pea", "wheat", "pea", "wheat")))
##' (Z_SIS_wheat <- mkZSIS(dat, colIDfocal="focal", colStand="stand",
##'                        colS="species", species="wheat",
##'                        weight=0.5))
##' (Z_SIS_pea <- mkZSIS(dat, colIDfocal="focal", colStand="stand",
##'                      colS="species", species="pea",
##'                      weight=0.25))
##' @noRd
mkZSIS <- function(df, colS, species, weight, colIDfocal = "focal", colStand = "stand",
                   levOneComp = "sole", levSevComps = "inter") {
  stopifnot(
    is.data.frame(df),
    all(sort(unique(as.character(df[[colStand]]))) %in%
      sort(c(levOneComp, levSevComps))),
    all(!is.na(df[[colStand]]))
  )

  idx <- which(df[[colS]] == species)
  genos <- sort(unique(as.character(df[[colIDfocal]][idx])))
  nbGenos <- length(genos)
  Z_SIS <- matrix(0, nrow = nrow(df), ncol = nbGenos)
  colnames(Z_SIS) <- genos

  for (geno in genos) {
    isPure <- (df[[colStand]] == levOneComp & df[[colIDfocal]] == geno)
    if (any(isPure)) {
      Z_SIS[which(isPure), geno] <- 1
    }
    isMixed <- (df[[colStand]] == levSevComps & df[[colIDfocal]] == geno)
    if (any(isMixed)) {
      Z_SIS[which(isMixed), geno] <- weight
    }
  }

  return(Z_SIS)
}

## ============================================================================
## code to fit models

##' Pseudo-R2
##'
##' Computes a pseudo-R2 as proposed by Efron (1978).
##' @param fit "merMod" object returned by \code{\link[lme4]{lmer}}
##' @return numeric
##' @author Timothee Flutre
##' @noRd
pseudoR2 <- function(fit) {
  stopifnot(inherits(fit, "lmerMod"))

  y <- fit@resp$y
  y_hat <- fitted(fit) # same as fit@resp$mu, but more generic
  y0 <- mean(y) # slightly different than fit@beta[1]

  residuals <- y - y_hat
  SSR <- sum(residuals^2)
  SST <- sum((y - y0)^2)
  1 - SSR / SST
}

##' Model summary
##'
##' Summarizes the GMA-SMA models.
##' @param fits named list of "merMod" objects returned by \code{\link[lme4]{lmer}}, as is done by \code{\link{fitGMASMA}} with \code{pkg="lme4"}
##' @return matrix
##' @seealso \code{\link{fitGMASMA}}
##' @author Timothee Flutre
##' @examples
##' ## see the first vignette
##' @export
summarizeGMASMAs <- function(fits) {
  mods <- paste0("mod", c("1", "2p", "2", "2pp", "3", "3p"))
  stopifnot(
    is.list(fits),
    !is.null(names(fits)),
    all(names(fits) %in% mods)
  )

  ## structure the output
  out <- list()
  colNs <- c("AIC", "R2", "RMSE", "var.GMA", "var.SMA", "var.SMAij", "var.SMAii", "var.err")
  for (mod in names(fits)) {
    out[[mod]] <- c()
    for (x in colNs) {
      out[[mod]][[x]] <- NA
    }
  }

  ## fill the output
  for (mod in names(fits)) {
    for (x in colNs) {
      if (x == "AIC") {
        out[[mod]][x] <- extractAIC(fits[[mod]])[2]
      } else if (x == "R2") {
        out[[mod]][x] <- pseudoR2(fits[[mod]])
      } else if (x == "RMSE") {
        out[[mod]][x] <- sqrt(mean(residuals(fits[[mod]])^2))
      } else if (x == "var.GMA") {
        tmp <- as.data.frame(lme4::VarCorr(fits[[mod]]))
        out[[mod]][x] <- tmp[which(tmp$grp == "GMA"), "vcov"]
      } else if (x == "var.SMA") {
        tmp <- as.data.frame(lme4::VarCorr(fits[[mod]]))
        if ("SMA" %in% tmp$grp) {
          out[[mod]][x] <- tmp[which(tmp$grp == "SMA"), "vcov"]
        }
      } else if (x == "var.SMAij") {
        tmp <- as.data.frame(lme4::VarCorr(fits[[mod]]))
        if ("SMA_ij" %in% tmp$grp) {
          out[[mod]][x] <- tmp[which(tmp$grp == "SMA_ij"), "vcov"]
        }
      } else if (x == "var.SMAii") {
        tmp <- as.data.frame(lme4::VarCorr(fits[[mod]]))
        if ("SMA_ii" %in% tmp$grp) {
          out[[mod]][x] <- tmp[which(tmp$grp == "SMA_ii"), "vcov"]
        }
      } else if (x == "var.err") {
        tmp <- as.data.frame(lme4::VarCorr(fits[[mod]]))
        out[[mod]][x] <- tmp[which(tmp$grp == "Residual"), "vcov"]
      }
    }
  }

  out <- do.call(rbind, out)

  return(out)
}

##' Fit GMA-SMA models
##'
##' Fits GMA-SMA models.
##' @param formFix formula whose right-hand side should only include the fixed effects, the random effects being deduced based on the names of \code{"listZ"}
##' @param data data frame; missing data will be removed
##' @param listZ named list of design matrices for the random effects; its names should be \code{"GMA"} (compulsory) and, optionally, one of \code{"SMA"}, \code{"SMA_ij"} or \code{"SMA_ii"}
##' @param pkg package used to fit the model among lme4, MM4LMM, TMB, INLA and breedR; lme4 and INLA do not use \code{"listVCov"}; breedR does not use \code{"contrasts"}
##' @param listVCov named list of variance-covariance matrices for the random effects; the names of this list should be the same as the names of \code{"listZ"}
##' @param contrasts  see \code{\link[stats:model.matrix]{stats::model.matrix()}}
##' @param REML logical
##' @param ... additional arguments specific to each package
##' @return depends on the chosen package
##' @seealso \code{\link{mkZGMA}}, \code{\link{mkZSMA}}
##' @author Timothee Flutre
##' @examples
##' ## see the first and second vignettes
##' @export
fitGMASMA <- function(formFix, data, listZ, pkg = "lme4", listVCov, contrasts = NULL,
                      REML = TRUE, ...) {
  stopifnot(pkg %in% c("lme4", "MM4LMM", "TMB", "INLA", "breedR"))
  if (pkg == "lme4") {
    lmerGMASMA(formFix, data, listZ, contrasts, REML)
  } else if (pkg == "MM4LMM") {
    mmGMASMA(formFix, data, listZ, listVCov, contrasts, REML)
  } else if (pkg == "TMB") {
    tmbGMASMA(formFix, data, listZ, listVCov, REML, contrasts, ...)
  } else if (pkg == "INLA") {
    inlaGMASMA(formFix, data, listZ, contrasts, ...)
  } else if (pkg == "breedR") {
    suppressMessages(bfGMASMA(formFix, data, listZ, listVCov))
  }
}

##' Fit a GMA-SMA model with lme4
##'
##' Fits a GMA-SMA model with \href{https://cran.r-project.org/package=lme4}{lme4}.
##' @param formFix see \code{\link[lme4:lmer]{lme4::lmer()}}
##' @param data see \code{\link[lme4:lmer]{lme4::lmer()}}; missing data will be removed
##' @param listZ named list of design matrices; its names should be \code{"GMA"} (compulsory) and, optionally, one of \code{"SMA"}, \code{"SMA_ij"} or \code{"SMA_ii"}
##' @param REML see \code{\link[lme4:lmer]{lme4::lmer()}}
##' @param contrasts see \code{\link[lme4:lmer]{lme4::lmer()}}
##' @return merMod object
##' @seealso \code{\link{mkZGMA}}, \code{\link{mkZSMA}}, \code{\link{fitGMASMA}}
##' @author Timothee Flutre
##' @examples
##' ## generate fake data
##' nbGenos <- 25
##' genos <- sprintf("g%02i", 1:nbGenos)
##' pairs <- t(combn(x=genos, m=2))
##' stands <- paste(pairs[,1], pairs[,2], sep="_")
##' nbBlocks <- 3
##' blocks <- LETTERS[1:nbBlocks]
##' dat <- do.call(rbind, lapply(blocks, function(block){
##'   cbind(stands=as.data.frame(stands, stringsAsFactors=TRUE),
##'         block=as.factor(block))
##' }))
##' listContr <- list(block="contr.sum")
##' X <- model.matrix(~ 1 + block, data=dat, contrasts=listContr)
##' Z_GMA <- mkZGMA(df=dat, col="stands", sep="_")
##' image(t(Z_GMA)[, nrow(Z_GMA):1], main="Z_GMA", axes=FALSE)
##' truth <- list("intercept"=100,
##'               "var_GMA"=10,
##'               "var_error"=1)
##' set.seed(1234)
##' truth[["blockEffs"]] <- sample(x=c(-1,1), size=nbBlocks - 1, replace=TRUE) *
##'   rnorm(n=nbBlocks - 1, mean=3, sd=5)
##' truth[["GMAs"]] <- rnorm(n=nbGenos, mean=0, sd=sqrt(truth$var_GMA))
##' truth[["errors"]] <- rnorm(n=nrow(dat), mean=0, sd=sqrt(truth$var_error))
##' y <- X %*% c(truth$intercept, truth$blockEffs) +
##'   Z_GMA %*% truth$GMAs +
##'   truth$errors
##' dat$pheno <- y[,1]
##' if(FALSE){
##'   hist(dat$pheno, las=1, main="Simulated data")
##'   boxplot(pheno ~ block, data=dat, las=1, main="Simulated data")
##' }
##'
##' ## fit the model
##' myformFix <- pheno ~ 1 + block
##' listZ <- list("GMA"=Z_GMA)
##' system.time(
##'   fit <- lmerGMASMA(myformFix, dat, listZ, listContr))
##'
##' ## check the results
##' BLUEs <- fixef(fit)
##' data.frame("true"=truth$intercept,
##'            "estim"=BLUEs["(Intercept)"])
##' data.frame("true"=truth$blockEffs,
##'            "estim"=BLUEs[grep("block", names(BLUEs))])
##' estV <- as.data.frame(lme4::VarCorr(fit))
##' data.frame("true"=c(truth$var_GMA, truth$var_error),
##'            "estim"=estV$vcov,
##'            row.names=c("GMA","error"))
##' BLUPs <- ranef(fit)
##' cor(truth$GMAs, BLUPs$GMA[,"(Intercept)"])
##' plot(BLUPs$GMA[,"(Intercept)"], truth$GMAs,
##'      xlab="BLUPs(GMA)", ylab="true GMAs",
##'      main="Accuracy with lme4", las=1, pch=19)
##' abline(a=0, b=1, v=0, h=0, lty=2)
##' abline(lm(truth$GMAs ~ BLUPs$GMA[,"(Intercept)"]), col="red")
##' @noRd
lmerGMASMA <- function(formFix, data, listZ, contrasts = NULL, REML = TRUE) {
  SMAtypes <- c("SMA", "SMA_ij", "SMA_ii")
  withSMAs <- setNames(rep(FALSE, 3), SMAtypes)

  ## checks
  stopifnot(
    requireNamespace("Matrix", quietly = TRUE),
    requireNamespace("lme4", quietly = TRUE),
    is.data.frame(data),
    is.list(listZ),
    !is.null(names(listZ)),
    "GMA" %in% names(listZ),
    is.matrix(listZ$GMA)
  )
  for (SMAtype in SMAtypes) {
    withSMAs[SMAtype] <- SMAtype %in% names(listZ)
    if (withSMAs[SMAtype]) {
      stopifnot(is.matrix(listZ[[SMAtype]]))
    }
  }

  ## discard missing data
  response <- as.character(formFix)[2]
  if (anyNA(data[[response]])) {
    idxNotNA <- which(!is.na(data[[response]]))
    msg <- paste0("removed missing data (", nrow(data) - length(idxNotNA), ")")
    warning(msg)
    data <- droplevels(data[idxNotNA, ])
    listZ <- lapply(listZ, function(Z) {
      Z[idxNotNA, ]
    })
  }

  ## update the input data frame
  if (!"GMA" %in% colnames(data)) {
    data <- cbind(data, "GMA" = rep(0, nrow(data)))
  }
  for (SMAtype in SMAtypes) {
    if (withSMAs[SMAtype]) {
      if (!SMAtype %in% colnames(data)) {
        data[[SMAtype]] <- data[["GMA"]]
      }
    }
  }

  ## update the formula
  formTxt <- as.character(formFix)
  formTxt <- paste0(
    formTxt[2], " ~ ", formTxt[3],
    " + (1|GMA)"
  )
  for (SMAtype in SMAtypes) {
    if (withSMAs[SMAtype]) {
      formTxt <- paste0(formTxt, " + (1|", SMAtype, ")")
    }
  }
  formula <- as.formula(formTxt)
  rm(formTxt)

  ## input design matrices
  myblist <- list()
  myblist[["GMA"]] <- list(
    "ff" = factor(colnames(listZ$GMA)),
    "sm" = Matrix::Matrix(t(listZ$GMA), sparse = TRUE),
    "nl" = as.integer(ncol(listZ$GMA)),
    "cnms" = "(Intercept)"
  )
  for (SMAtype in SMAtypes) {
    if (withSMAs[SMAtype]) {
      myblist[[SMAtype]] <- list(
        "ff" = factor(colnames(listZ[[SMAtype]])),
        "sm" = Matrix::Matrix(t(listZ[[SMAtype]]), sparse = TRUE),
        "nl" = as.integer(ncol(listZ[[SMAtype]])),
        "cnms" = "(Intercept)"
      )
    }
  }
  stopifnot(all(names(myblist) ==
    barnames(lme4::findbars(RHSForm(formula)))))

  ## model fit
  fit <- lmerZ(
    formula = formula, data = data, REML = REML, myblist = myblist,
    na.action = stats::na.fail, contrasts = contrasts,
    control = lme4::lmerControl(
      check.nobs.vs.nlev = "ignore",
      check.nobs.vs.nRE = "ignore"
    )
  )
  ## TODO: renvoie un warning dans mkLmerDevfun -> ... -> ave(), mais a priori pas grave

  return(fit)
}

##' Fit all three GMA-SMA models from Forst et al (2019) with lme4
##'
##' Fit all three GMA-SMA models from \href{https://doi.org/10.1016/j.fcr.2019.107571}{Forst et al (2019)} with \href{https://cran.r-project.org/package=lme4}{lme4}.
##' @param formFix see \code{\link[lme4:lmer]{lme4::lmer()}} with neither \code{"(1|GMA)"} nor \code{"(1|SMA)"}
##' @param data see \code{\link[lme4:lmer]{lme4::lmer()}}; missing data will be removed
##' @param listZ named list of three design matrices; its names should be \code{"GMA"}, \code{"SMA2"} and \code{"SMA3"}, all being compulsory
##' @param REML see \code{\link[lme4:lmer]{lme4::lmer()}}
##' @param contrasts see \code{\link[lme4:lmer]{lme4::lmer()}}
##' @return list of two components, the first with merMod objects and the second with a summary
##' @seealso \code{\link{mkZGMA}}, \code{\link{mkZSMA}}, \code{\link{fitGMASMA}}
##' @author Timothee Flutre
##' @examples
##' ## see the example in ?lmerGMASMA
##' @noRd
lmerGMASMAs <- function(formFix, data, listZ, contrasts = NULL, REML = TRUE) {
  stopifnot(
    is.list(listZ),
    length(listZ) >= 3,
    !is.null(names(listZ)),
    all(c("GMA", "SMA2", "SMA3") %in% names(listZ))
  )

  out <- list()

  ## make formula of each model
  formChar <- as.character(formFix)
  response <- formChar[2]
  RHS <- formChar[3] # right-hand side
  formulaGMA <- stats::as.formula(paste(response, "~", RHS, "+", "(1|GMA)"))
  formulaGMASMA <- stats::as.formula(paste(
    response, "~", RHS,
    "+", "(1|GMA)", "+", "(1|SMA)"
  ))

  ## make listZ of each model
  listZGMA <- list("GMA" = listZ$GMA)
  listZGMASMA2 <- list("GMA" = listZ$GMA, "SMA" = listZ$SMA2)
  listZGMASMA3 <- list("GMA" = listZ$GMA, "SMA" = listZ$SMA3)

  ## fit each model
  mods <- list()
  mods[["mod1"]] <- lmerGMASMA(
    formFix, data, listZGMA,
    contrasts, REML
  )
  mods[["mod2"]] <- lmerGMASMA(
    formFix, data, listZGMASMA2,
    contrasts, REML
  )
  mods[["mod3"]] <- lmerGMASMA(
    formFix, data, listZGMASMA3,
    contrasts, REML
  )
  out$objs <- mods

  ## summarize each model
  sry <- c()
  for (x in c("AIC", "LRT.pval", "R2", "var.GMA", "var.SMA", "var.err")) {
    for (mod in paste0("mod", 1:3)) {
      sry[[paste0(x, ".", mod)]] <- NA
    }
  }
  for (x in c("AIC", "LRT.pval", "R2", "var.GMA", "var.SMA", "var.err")) {
    for (mod in paste0("mod", 1:3)) {
      name <- paste0(x, ".", mod)
      if (x == "AIC") {
        sry[name] <- stats::extractAIC(mods[[mod]])[2]
      } else if (x == "LRT.pval") {
        if (mod == "mod2") {
          tmp <- suppressMessages(
            stats::anova(mods[["mod1"]], mods[["mod2"]])
          )
          sry[name] <- tmp[2, "Pr(>Chisq)"]
        } else if (mod == "mod3") {
          tmp <- suppressMessages(
            stats::anova(mods[["mod1"]], mods[["mod3"]])
          )
          sry[name] <- tmp[2, "Pr(>Chisq)"]
        }
      } else if (x == "R2") {
        sry[name] <- myRsq(mods[[mod]])
      } else if (x == "var.GMA") {
        tmp <- as.data.frame(lme4::VarCorr(mods[[mod]]))
        sry[name] <- tmp[which(tmp$grp == "GMA"), "vcov"]
      } else if (x == "var.SMA") {
        tmp <- as.data.frame(lme4::VarCorr(mods[[mod]]))
        if ("SMA" %in% tmp$grp) {
          sry[name] <- tmp[which(tmp$grp == "SMA"), "vcov"]
        }
      } else if (x == "var.err") {
        tmp <- as.data.frame(lme4::VarCorr(mods[[mod]]))
        sry[name] <- tmp[which(tmp$grp == "Residual"), "vcov"]
      }
    }
  }
  out[["sry"]] <- sry

  return(out)
}

##' Fit a GMA-SMA model with MM4LMM
##'
##' Fits a GMA-SMA model with \href{https://cran.r-project.org/package=MM4LMM}{MM4LMM}.
##' @param formFix formula whose right-hand side should only include the fixed effects, the random effects being deduced based on the names of \code{"listZ"}
##' @param data data data frame; missing data will be removed
##' @param listZ named list of design matrices for the random effects; its names should be \code{"GMA"} (compulsory) and, optionally, one of \code{"SMA"}, \code{"SMA_ij"} or \code{"SMA_ii"}; the matrix for the errors will be added automatically
##' @param listVCov named list of variance-covariance matrices for the random effects; the names of this list should be the same as the names of \code{"listZ"}; the matrix for the errors will be added automatically
##' @param contrasts see \code{\link[stats:model.matrix]{stats::model.matrix()}}
##' @param REML TRUE or FALSE
##' @return list
##' @seealso \code{\link{mkZGMA}}, \code{\link{mkZSMA}}, \code{\link{fitGMASMA}}
##' @author Timothee Flutre
##' @examples
##' \dontrun{
##' ## generate fake data
##' nbGenos <- 25
##' genos <- sprintf("g%02i", 1:nbGenos)
##' pairs <- t(combn(x=genos, m=2))
##' stands <- paste(pairs[,1], pairs[,2], sep="_")
##' nbBlocks <- 3
##' blocks <- LETTERS[1:nbBlocks]
##' dat <- do.call(rbind, lapply(blocks, function(block){
##'   cbind(stands=as.data.frame(stands, stringsAsFactors=TRUE),
##'         block=as.factor(block))
##' }))
##' listContr <- list(block="contr.sum")
##' X <- model.matrix(~ 1 + block, data=dat, contrasts=listContr)
##' Z_GMA <- mkZGMA(df=dat, col="stands", sep="_")
##' truth <- list("intercept"=100,
##'               "var_GMA"=10,
##'               "var_error"=1)
##' set.seed(1234)
##' truth[["blockEffs"]] <- sample(x=c(-1,1), size=nbBlocks - 1, replace=TRUE) *
##'   rnorm(n=nbBlocks - 1, mean=3, sd=5)
##' truth[["GMAs"]] <- rnorm(n=nbGenos, mean=0, sd=sqrt(truth$var_GMA))
##' truth[["errors"]] <- rnorm(n=nrow(dat), mean=0, sd=sqrt(truth$var_error))
##' y <- X %*% c(truth$intercept, truth$blockEffs) +
##'   Z_GMA %*% truth$GMAs +
##'   truth$errors
##' dat$pheno <- y[,1]
##' if(FALSE){
##'   hist(dat$pheno, las=1, main="Simulated data")
##'   boxplot(pheno ~ block, data=dat, las=1, main="Simulated data")
##' }
##'
##' ## fit the model
##' myformFix <- pheno ~ 1 + block
##' listZ <- list("GMA"=Z_GMA, "Error"=diag(nrow(dat)))
##' listVCov <- list("GMA"=diag(ncol(Z_GMA)), "Error"=diag(nrow(dat)))
##' system.time(
##'   fit <- mmGMASMA(myformFix, dat, listZ, listVCov, listContr))
##'
##' ## check the results
##' BLUEs <- fit$Beta
##' data.frame("true"=truth$intercept,
##'            "estim"=BLUEs["(Intercept)"])
##' data.frame("true"=truth$blockEffs,
##'            "estim"=BLUEs[grep("block", names(BLUEs))])
##' estV <- fit$Sigma2
##' data.frame("true"=c(truth$var_GMA, truth$var_error),
##'            "estim"=estV,
##'            row.names=c("GMA","error"))
##' BLUPs <- fit$BLUPs
##' cor(truth$GMAs, BLUPs$GMA[,1])
##' plot(BLUPs$GMA[,1], truth$GMAs,
##'      xlab="BLUPs(GMA)", ylab="true GMAs",
##'      main="Accuracy with MM4LMM", las=1, pch=19)
##' abline(a=0, b=1, v=0, h=0, lty=2)
##' abline(lm(truth$GMAs ~ BLUPs$GMA[,1]), col="red")
##' }
##' @noRd
mmGMASMA <- function(formFix, data, listZ, listVCov, contrasts = NULL, REML = TRUE) {
  stopifnot(
    requireNamespace("MM4LMM", quietly = TRUE),
    is.list(listZ),
    !is.null(names(listZ)),
    "GMA" %in% names(listZ),
    is.list(listVCov),
    !is.null(names(listVCov)),
    all(names(listVCov) == names(listZ))
  )

  out <- NULL

  ## input matrices:
  ## * design matrix of fixed effects:
  formChar <- as.character(formFix)
  response <- formChar[2]
  RHS <- formChar[3] # right-hand side
  X_C <- model.matrix(as.formula(paste("~", RHS)),
    data = data, contrasts.arg = contrasts
  )
  ## * design matrix of errors:
  listZ$Error <- diag(nrow(data))
  ## * vcov matrix of errors:
  listVCov$Error <- diag(nrow(data))

  ## model fit
  fit <- MM4LMM::MMEst(
    Y = data[[response]],
    Cofactor = X_C,
    ZList = listZ,
    VarList = listVCov,
    Method = ifelse(REML, "Reml", "ML")
  )
  if (FALSE) { # debug with example above:
    out <- MM4LMM::MMEst(Y = dat$pheno, Cofactor = model.matrix(~ 1 + block, dat), VarList = list("GMA" = diag(ncol(Z_GMA)), "Error" = diag(1, nrow(dat))), ZList = list("GMA" = Z_GMA, "Error" = diag(1, nrow(dat))))
  }
  out <- fit$NullModel

  ## BLUP computations
  tmp <- MM4LMM::MMBlup(
    Y = data[[response]],
    Cofactor = X_C,
    ZList = listZ,
    VarList = listVCov,
    ResMM = fit
  )
  for (i in 1:length(tmp)) {
    rndEff <- names(tmp)[[i]]
    if (is.null(rownames(tmp[[rndEff]]))) {
      rownames(tmp[[rndEff]]) <- colnames(listZ[[rndEff]])
    }
  }
  out$BLUPs <- tmp

  return(out)
}

##' Fit all three GMA-SMA models with MM4LMM
##'
##' Fits all three GMA-SMA models of \href{https://doi.org/10.1016/j.fcr.2019.107571}{Forst et al (2019)} with \href{https://cran.r-project.org/package=MM4LMM}{MM4LMM}.
##' @param formFix see \code{\link[lme4:lmer]{lme4::lmer()}} with neither \code{"(1|GMA)"} nor \code{"(1|SMA)"}
##' @param data see \code{\link[lme4:lmer]{lme4::lmer()}}; missing data will be removed
##' @param listZ see argument \code{"ZList"} of \code{\link[MM4LMM:MMEst]{MM4LMM::MMEst()}}; the names of this list should contain at least \code{"GMA"} and \code{"Error"} and, optionally, \code{"SMA2"} and \code{"SMA3"}
##' @param listVCov see argument \code{"VarList"} of \code{\link[MM4LMM:MMEst]{MM4LMM::MMEst()}}; the names of this list should be the same as the names of \code{"listZ"}
##' @param contrasts see \code{\link[stats:model.matrix]{stats::model.matrix()}}
##' @param REML see \code{\link[lme4:lmer]{lme4::lmer()}}
##' @param nbCores see argument \code{"mc.cores"} of \code{\link[parallel:mclapply]{paralell::mclapply()}}
##' @return list of two components, the first with list objects and the second with a summary
##' @seealso \code{\link{mkZGMA}}, \code{\link{mkZSMA}}, \code{\link{mmGMASMA}}
##' @author Timothee Flutre
##' @examples
##' ## see the example in ?mmGMASMA
##' @noRd
mmGMASMAs <- function(formFix, data, listZ, listVCov, contrasts = NULL,
                      REML = TRUE, nbCores = 1) {
  stopifnot(
    is.list(listZ),
    length(listZ) >= 3,
    !is.null(names(listZ)),
    all(c("GMA", "SMA2", "SMA3", "Error") %in% names(listZ)),
    is.list(listVCov),
    all(names(listVCov) == names(listZ))
  )

  out <- list()

  ## fit each model
  mods <- parallel::mclapply(1:3, function(i) {
    if (i == 1) {
      mmGMASMA(formFix, data,
        listZ[c("GMA", "Error")],
        listVCov[c("GMA", "Error")],
        contrasts = contrasts, REML = REML
      )
    } else if (i == 2) {
      mmGMASMA(formFix, data,
        list("GMA" = listZ$GMA, "SMA" = listZ$SMA2, "Error" = listZ$Error),
        list("GMA" = listVCov$GMA, "SMA" = listVCov$SMA2, "Error" = listVCov$Error),
        contrasts = contrasts, REML = REML
      )
    } else if (i == 3) {
      mmGMASMA(formFix, data,
        list("GMA" = listZ$GMA, "SMA" = listZ$SMA3, "Error" = listZ$Error),
        list("GMA" = listVCov$GMA, "SMA" = listVCov$SMA3, "Error" = listVCov$Error),
        contrasts = contrasts, REML = REML
      )
    }
  }, mc.cores = nbCores)
  names(mods) <- paste0("mod", 1:3)
  out$objs <- mods

  ## TODO: summarize each model
  sry <- c()
  for (x in c("AIC", "LRT.pval", "R2", "var.GMA", "var.SMA", "var.err")) {
    for (mod in paste0("mod", 1:3)) {
      sry[[paste0(x, ".", mod)]] <- NA
    }
  }

  return(out)
}

##' Fit a DGE-IGE model with MM4LMM
##'
##' Fits a DGE-IGE model with \href{https://cran.r-project.org/package=MM4LMM}{MM4LMM}.
##' @param formula see \code{\link[stats:lm]{stats::lm()}} only with fixed effects
##' @param data see \code{\link[stats:lm]{stats::lm()}}
##' @param listZ see argument \code{"ZList"} of \code{\link[MM4LMM:MMEst]{MM4LMM::MMEst()}}; the names of this list should contain \code{"DGE"}, \code{"IGE"} and \code{"Error"}
##' @param listVCov see argument \code{"VarList"} of \code{\link[MM4LMM:MMEst]{MM4LMM::MMEst()}}; the names of this list should be the same as the names of \code{"listZ"}
##' @param REML if TRUE, the restricted log-likelihood is maximized, otherwise it is the log-likelihood
##' @param contrasts.arg see \code{\link[stats:model.matrix]{stats::model.matrix()}}
##' @param CritVar see \code{\link[MM4LMM:MMEst]{MM4LMM::MMEst()}}
##' @param CritLogLik see \code{\link[MM4LMM:MMEst]{MM4LMM::MMEst()}}
##' @param MaxIter see \code{\link[MM4LMM:MMEst]{MM4LMM::MMEst()}}
##' @param verbose verbosity level
##' @return list
##' @author Tristan Mary-Huard (aut), Timothee Flutre (ctb)
##' @examples
##' \dontrun{
##' ## generate fake data
##' nbGenos <- 25
##' genos <- sprintf("g%02i", 1:nbGenos)
##' pairs <- t(combn(x=genos, m=2))
##' nbMixtures <- nrow(pairs)
##' stands <- paste(pairs[,1], pairs[,2], sep="_")
##' nbBlocks <- 3
##' blocks <- LETTERS[1:nbBlocks]
##' dat <- data.frame(stand=rep(stands, each=2),
##'                   block=rep(blocks, each=nbMixtures * 2),
##'                   DGE=NA,
##'                   IGE=NA,
##'                   stringsAsFactors=TRUE)
##' idx1 <- seq(1, nrow(dat)-1, by=2)
##' dat$DGE[idx1] <- sapply(strsplit(as.character(dat$stand)[idx1], "_"), `[`, 1)
##' dat$IGE[idx1] <- sapply(strsplit(as.character(dat$stand)[idx1], "_"), `[`, 2)
##' idx2 <- seq(2, nrow(dat), by=2)
##' dat$DGE[idx2] <- sapply(strsplit(as.character(dat$stand)[idx2], "_"), `[`, 2)
##' dat$IGE[idx2] <- sapply(strsplit(as.character(dat$stand)[idx2], "_"), `[`, 1)
##' dat$DGE <- factor(dat$DGE, levels=genos)
##' dat$IGE <- factor(dat$IGE, levels=genos)
##' listContr <- list(block="contr.sum")
##' X <- model.matrix(~ 1 + block, data=dat, contrasts=listContr)
##' tmp <- mkZintraspe(dat, colIDfocal="DGE", colIDneighbors="IGE"); Z_D <- tmp$Z_DBV; Z_I <- tmp$Z_SBV
##' truth <- list("intercept"=100,
##'               "var_DGE"=10,
##'               "var_IGE"=2,
##'               "cor_DGE-IGE"=-0.8,
##'               "var_error"=1)
##' set.seed(1234)
##' truth[["blockEffs"]] <- sample(x=c(-1,1), size=nbBlocks - 1, replace=TRUE) *
##'   rnorm(n=nbBlocks - 1, mean=3, sd=5)
##' Sigma <- matrix(c(truth$var_DGE, NA, NA, truth$var_IGE), 2, 2)
##' Sigma[1,2] <- Sigma[2,1] <- truth[["cor_DGE-IGE"]] * sqrt(Sigma[1,1] * Sigma[2,2])
##' tmp <- mvtnorm::rmvnorm(n=nbGenos, sigma=Sigma)
##' truth[["DGEs"]] <- tmp[,1]
##' truth[["IGEs"]] <- tmp[,2]
##' plot(truth[["DGEs"]], truth[["IGEs"]], las=1, main="Simulated genetic effects")
##' abline(h=0, v=0, lty=2)
##' abline(lm(truth[["IGEs"]] ~ truth[["DGEs"]]), col="red")
##' truth[["errors"]] <- rnorm(n=nrow(dat), mean=0, sd=sqrt(truth$var_error))
##' y <- X %*% c(truth$intercept, truth$blockEffs) +
##'   Z_D %*% truth[["DGEs"]] +
##'   Z_I %*% truth[["IGEs"]] +
##'   truth$errors
##' dat$pheno <- y[,1]
##' head(dat)
##' if(FALSE){
##'   hist(dat$pheno, las=1, main="Simulated data")
##'   boxplot(pheno ~ block, data=dat, las=1, main="Simulated data")
##' }
##'
##' ## fit the model
##' myformula <- pheno ~ 1 + block
##' listZ <- list("DGE"=Z_D, "IGE"=Z_I, "Error"=diag(nrow(dat)))
##' listVCov <- list("DGE"=diag(ncol(Z_D)), "IGE"=diag(ncol(Z_I)), "Error"=diag(nrow(dat)))
##' system.time(
##'   fit <- mmDGEIGE(myformula, dat, listZ, listVCov, contrasts.arg=listContr))
##'
##' ## check the results
##' BLUEs <- fit$Beta
##' data.frame("true"=c(truth$intercept, truth$blockEffs),
##'            "estim"=BLUEs,
##'            row.names=names(BLUEs))
##' estV <- fit$Sigma2
##' estV <- c(estV, "CorDI"=as.numeric(estV["CovDI"] / sqrt(estV["DGE"] * estV["IGE"])))
##' data.frame("true"=c(truth$var_DGE, truth$var_IGE, truth$`cor_DGE-IGE`, truth$var_error),
##'            "estim"=estV[c("DGE", "IGE", "CorDI", "Error")],
##'            row.names=c("DGE","IGE","CorDI","error"))
##' BLUPs <- fit$BLUPs
##' cor(truth$DGEs, BLUPs$DGE[,1])
##' cor(truth$IGEs, BLUPs$IGE[,1])
##' plot(BLUPs$DGE[,1], truth$DGEs,
##'      xlab="BLUPs(DGE)", ylab="true DGEs",
##'      main="Accuracy with MM4LMM", las=1, pch=19)
##' abline(a=0, b=1, v=0, h=0, lty=2)
##' abline(lm(truth$DGEs ~ BLUPs$DGE[,1]), col="red")
##' plot(BLUPs$IGE[,1], truth$IGEs,
##'      xlab="BLUPs(IGE)", ylab="true IGEs",
##'      main="Accuracy with MM4LMM", las=1, pch=19)
##' abline(a=0, b=1, v=0, h=0, lty=2)
##' abline(lm(truth$IGEs ~ BLUPs$IGE[,1]), col="red")
##' }
##' @noRd
mmDGEIGE <- function(formula, data, listZ, listVCov, REML = TRUE,
                     contrasts.arg = NULL,
                     CritVar = 0.001, CritLogLik = 0.001, MaxIter = 100,
                     verbose = FALSE) {
  ## to debug
  ## formula <- myformula; data <- dat; REML=TRUE; contrasts.arg <- listContr; CritVar=0.001; CritLogLik=0.001; MaxIter=100
  stopifnot(
    requireNamespace("MM4LMM", quietly = TRUE),
    is.list(listZ),
    all(c("DGE", "IGE", "Error") %in% names(listZ)),
    all(names(listVCov) == names(listZ))
  )

  out <- list()

  formChar <- as.character(formula)
  response <- formChar[2]
  RHS <- formChar[3] # right-hand side
  y <- data[[response]]
  X_C <- model.matrix(as.formula(paste("~", RHS)),
    data = data, contrasts.arg = contrasts.arg
  )
  listZ <- lapply(listZ, as.matrix)
  Z_D <- listZ[["DGE"]]
  Z_I <- listZ[["IGE"]]
  listVCov <- lapply(listVCov, as.matrix)

  if (verbose) {
    msg <- "initialization: fit the model without cov(DGE,IGE)"
    print(msg)
  }
  initFit <- MM4LMM::MMEst(
    Y = data[[response]],
    Cofactor = X_C,
    ZList = listZ,
    VarList = listVCov,
    Method = ifelse(REML, "Reml", "ML"),
    Henderson = FALSE, CritVar = CritVar
  )
  initBLUPs <- MM4LMM::MMBlup(
    Y = data[[response]],
    Cofactor = X_C,
    ZList = listZ,
    VarList = listVCov,
    ResMM = initFit
  )
  initVars <- initFit$NullModel$Sigma2
  initCorDI <- as.numeric(cor(initBLUPs[["DGE"]], initBLUPs[["IGE"]]))
  initCovDI <- initCorDI * sqrt(initVars[["DGE"]] * initVars[["IGE"]])
  signCovDI <- sign(initCovDI)
  if (verbose) {
    print(initVars)
    print(paste0("cor(BLUPs(DGE),BLUPs(IGE)) = ", round(initCorDI, 2)))
  }

  if (verbose) {
    msg <- "loop: re-fit the model with cov(DGE,IGE)"
    print(msg)
  }
  initVars <- c(
    DGE = initVars[["DGE"]] - initCovDI, # var_DGE
    IGE = initVars[["IGE"]] - initCovDI, # var_IGE
    CovDI = abs(initCovDI), # cov_DGE-IGE
    Error = initVars[["Error"]]
  ) # var_Err
  listVCov <- list(
    DGE = listVCov[["DGE"]],
    IGE = listVCov[["IGE"]],
    CovDI = listVCov[["DGE"]],
    Error = listVCov[["Error"]]
  )
  MinVar <- which.min(c(initVars["DGE"], initVars["IGE"]))
  Alpha <- 1
  NullMinVar <- TRUE
  First <- TRUE
  iter <- 0
  loglik <- NA
  while (NullMinVar) {
    if (verbose) {
      print(paste0("start iter ", iter + 1))
    }
    Coef1 <- (MinVar == 1) * signCovDI * Alpha + (MinVar == 2) * 1
    Coef2 <- (MinVar == 2) * signCovDI * Alpha + (MinVar == 1) * 1
    listZ <- list(
      DGE = Z_D,
      IGE = Z_I,
      CovDI = Coef1 * Z_D + Coef2 * Z_I,
      Error = listZ[["Error"]]
    )
    fit <- MM4LMM::MMEst(
      Y = data[[response]],
      Cofactor = X_C,
      ZList = listZ,
      VarList = listVCov,
      Henderson = FALSE,
      CritVar = CritVar,
      Init = initVars
    )
    if (verbose) {
      cat("loglik:", fit$NullModel$`LogLik (Reml)`, "\n")
    }
    if (!is.na(loglik)) {
      if (abs(fit$NullModel$`LogLik (Reml)` - loglik) < CritLogLik) {
        if (verbose)
          cat("break because of CritLogLik at iter ", iter)
        break
      }
    }
    loglik <- fit$NullModel$`LogLik (Reml)`
    fit.var <- fit$NullModel$Sigma2
    if (fit.var[MinVar] < 1e-4) {
      Alpha <- 0.9 * Alpha
    } else {
      if (fit.var[1:2][-MinVar] < 1e-4) {
        if (First) {
          if (verbose)
            cat("somewhere in between ", Alpha, " and ", Alpha / 0.9, "...\n")
          First <- FALSE
        }
        Alpha <- 1.02 * Alpha
      } else {
        NullMinVar <- FALSE
      }
    }
    iter <- iter + 1
    if (iter > MaxIter) {
      if (verbose)
        cat("break because of MaxIter at iter ", iter)
      break
    }
  }
  if (verbose) {
    print(paste0("alpha = ", round(Alpha, 2)))
  }

  ## output
  out <- fit$NullModel
  out$BLUPs <- MM4LMM::MMBlup(
    Y = data[[response]],
    Cofactor = X_C,
    ZList = listZ,
    VarList = listVCov,
    ResMM = fit
  )
  out$Sigma2 <- c(
    "DGE" = as.numeric(fit$NullModel$Sigma2["DGE"] +
      ((MinVar == 1) * (Alpha**2) +
        (MinVar == 2)) * fit$NullModel$Sigma2["CovDI"]),
    "IGE" = as.numeric(fit$NullModel$Sigma2["IGE"] +
      ((MinVar == 2) * (Alpha**2) +
        (MinVar == 1)) * fit$NullModel$Sigma2["CovDI"]),
    "CovDI" = as.numeric(signCovDI * Alpha * fit$NullModel$Sigma2["CovDI"]),
    "Error" = as.numeric(fit$NullModel$Sigma2["Error"])
  )
  if (verbose) {
    print(out$Sigma2)
    corDI <- out$Sigma2["CovDI"] / sqrt(out$Sigma2["DGE"] * out$Sigma2["IGE"])
    print(paste0("cor(DGE,IGE) = ", round(corDI, 2)))
  }

  return(out)
}

##' Fit a GMA-SMA model with INLA
##'
##' Fits a GMA-SMA model with \href{https://www.r-inla.org/}{INLA}.
##' @param formFix see argument \code{"formula"} of \code{\link[INLA:inla]{INLA::inla()}}; ONLY include terms modeled as "fixed"
##' @param data see argument \code{"data"} of \code{\link[INLA:inla]{INLA::inla()}}
##' @param listZ named list of design matrices; the names of this list should contain at least \code{"GMA"} and, optionally, one of \code{"SMA"}, \code{"SMA_ij"} or \code{"SMA_ii"}
##' @param contrasts see \code{\link[INLA:inla]{INLA::inla()}}
##' @param priorPrecGma list specifying the prior on the GMA precision, e.g., a penalized complexity one (more information with \code{inla.doc("pc.prec")})
##' @param priorPrecErr list specifying the prior on the error precision, e.g., a log-Gamma (more information with \code{inla.doc("loggamma")})
##' @param priorPrecSma list specifying the prior on the SMA precision, e.g., a penalized complexity one (more information with \code{inla.doc("pc.prec")})
##' @return return value of \code{\link[INLA:inla]{INLA::inla()}}
##' @seealso \code{\link{mkZGMA}}, \code{\link{mkZSMA}}, \code{\link{fitGMASMA}}
##' @author Timothee Flutre
##' @examples
##' \dontrun{
##' ## generate fake data
##' nbGenos <- 25
##' genos <- sprintf("g%02i", 1:nbGenos)
##' pairs <- t(combn(x=genos, m=2))
##' stands <- paste(pairs[,1], pairs[,2], sep="_")
##' nbBlocks <- 3
##' blocks <- LETTERS[1:nbBlocks]
##' dat <- do.call(rbind, lapply(blocks, function(block){
##'   cbind(stands=as.data.frame(stands, stringsAsFactors=TRUE),
##'         block=as.factor(block))
##' }))
##' listContr <- list(block="contr.sum")
##' X <- model.matrix(~ 1 + block, data=dat, contrasts=listContr)
##' Z_GMA <- mkZGMA(df=dat, col="stands", sep="_")
##' truth <- list("intercept"=100,
##'               "var_GMA"=10,
##'               "var_error"=1)
##' set.seed(1234)
##' truth[["blockEffs"]] <- sample(x=c(-1,1), size=nbBlocks - 1, replace=TRUE) *
##'   rnorm(n=nbBlocks - 1, mean=3, sd=5)
##' truth[["GMAs"]] <- rnorm(n=nbGenos, mean=0, sd=sqrt(truth$var_GMA))
##' truth[["errors"]] <- rnorm(n=nrow(dat), mean=0, sd=sqrt(truth$var_error))
##' y <- X %*% c(truth$intercept, truth$blockEffs) +
##'   Z_GMA %*% truth$GMAs +
##'   truth$errors
##' dat$pheno <- y[,1]
##' if(FALSE){
##'   hist(dat$pheno, las=1, main="Simulated data")
##'   boxplot(pheno ~ block, data=dat, las=1, main="Simulated data")
##' }
##'
##' ## fit the model
##' myformFix <- pheno ~ 1 + block
##' listZ <- list("GMA"=Z_GMA)
##' fit <- inlaGMASMA(myformFix, dat, listZ, listContr)
##'
##' ## check the results
##' suppressPackageStartupMessages(library(INLA))
##' fit$summary.fixed
##' cbind(c(truth$intercept, truth$blockEffs), fit$summary.fixed[,"mean"])
##' fit$summary.hyperpar
##' names(fit$internal.marginals.hyperpar)
##' m <- fit$internal.marginals.hyperpar[[1]]
##' m.var <- inla.tmarginal(function(x) 1/exp(x), m)
##' inla.zmarginal(m.var) # mean, stdev, quantiles
##' inla.mmarginal(m.var) # mode
##' plot(inla.smarginal(m.var), type="l", main="var_error")
##' abline(v=truth$var_error, col="red")
##' m <- fit$internal.marginals.hyperpar[[2]]
##' m.var <- inla.tmarginal(function(x) 1/exp(x), m)
##' inla.zmarginal(m.var) # mean, stdev, quantiles
##' inla.mmarginal(m.var) # mode
##' plot(inla.smarginal(m.var), type="l", main="var_GMA")
##' abline(v=truth$var_GMA, col="red")
##' idx <- (nrow(dat)+1):nrow(fit$summary.random$ID_GMA)
##' cor(truth$GMAs, fit$summary.random$ID_GMA$mean[idx])
##' plot(fit$summary.random$ID_GMA$mean[idx], truth$GMAs, xlab="BLUPs(GMA)", ylab="true GMAs",
##'      main="Accuracy with INLA", las=1, pch=19)
##' abline(a=0, b=1, v=0, h=0, lty=2)
##' abline(lm(truth$GMAs ~ fit$summary.random$ID_GMA$mean[idx]), col="red")
##' }
##' @noRd
inlaGMASMA <- function(formFix, data, listZ, contrasts = NULL,
                       priorPrecGma = list(
                         prior = "pc.prec",
                         params = c("u" = 1, "alpha" = 0.01)
                       ),
                       priorPrecErr = list(
                         prior = "loggamma",
                         params = c("shape" = 1, "rate" = 1)
                       ),
                       priorPrecSma = list(
                         prior = "pc.prec",
                         params = c("u" = 1, "alpha" = 0.01)
                       )) {
  stopifnot(
    requireNamespace("INLA", quietly = TRUE),
    is.list(listZ),
    !is.null(names(listZ)),
    "GMA" %in% names(listZ)
  )

  data$ID_GMA <- 1:nrow(data)

  ## make the whole formula with the GMAs
  formChar <- as.character(formFix)
  response <- formChar[2]
  RHS <- formChar[3] # right-hand side
  form <- paste0(
    response, " ~ ", RHS,
    " + f(ID_GMA, model=\"z\", Z=listZ$GMA",
    ", hyper=list(prec=list(prior=\"", priorPrecGma$prior, "\"",
    ", param = c(", priorPrecGma$params[1]
  )
  for (i in seq_along(priorPrecGma$params)[-1]) {
    form <- paste0(form, ", ", priorPrecGma$params[i])
  }
  form <- paste0(form, "))))")

  ## if requested, add the SMAs
  isSMA <- grepl("SMA", names(listZ))
  if (any(isSMA)) {
    for (SMAtype in names(listZ)[which(isSMA)]) {
      ID_SMA <- paste0("ID_", SMAtype)
      data[[ID_SMA]] <- 1:nrow(data)
      form <- paste0(
        form,
        " + f(", ID_SMA, ", model=\"z\", Z=listZ$", SMAtype,
        ", hyper=list(prec=list(prior=\"", priorPrecSma$prior, "\"",
        ", param = c(", priorPrecSma$params[1]
      )
      for (i in seq_along(priorPrecSma$params)[-1]) {
        form <- paste0(form, ", ", priorPrecSma$params[i])
      }
      form <- paste0(form, "))))")
    }
  }

  ## fit the model
  form <- as.formula(form)
  fit <- INLA::inla(
    formula = form, data = data, contrasts = contrasts,
    control.family = list(hyper = list(prec = priorPrecErr)),
    control.compute = list(dic = TRUE)
  )

  return(fit)
}

##' Fit a DGE-IGE model with INLA
##'
##' Fits a DGE-IGE model with \href{https://www.r-inla.org/}{INLA}.
##' @param formula see argument \code{"formula"} of \code{\link[INLA:inla]{INLA::inla()}}
##' @param data see argument \code{"data"} of \code{\link[INLA:inla]{INLA::inla()}}
##' @param listZ named list of design matrices; the names of this list should contain at least \code{"DGE"} and \code{"IGE"}
##' @param contrasts see \code{\link[INLA:inla]{INLA::inla()}}
##' @param verbose verbosity level
##' @return return value of \code{\link[INLA:inla]{INLA::inla()}}
##' @seealso \code{\link{mkZintraspe}}, \code{\link{mmDGEIGE}}
##' @author Timothee Flutre
##' @examples
##' \dontrun{
##' ## generate fake data
##' nbGenos <- 25
##' genos <- sprintf("g%02i", 1:nbGenos)
##' pairs <- t(combn(x=genos, m=2))
##' nbMixtures <- nrow(pairs)
##' stands <- paste(pairs[,1], pairs[,2], sep="_")
##' nbBlocks <- 3
##' blocks <- LETTERS[1:nbBlocks]
##' dat <- data.frame(stand=rep(stands, each=2),
##'                   block=rep(blocks, each=nbMixtures * 2),
##'                   DGE=NA,
##'                   IGE=NA,
##'                   stringsAsFactors=TRUE)
##' idx1 <- seq(1, nrow(dat)-1, by=2)
##' dat$DGE[idx1] <- sapply(strsplit(as.character(dat$stand)[idx1], "_"), `[`, 1)
##' dat$IGE[idx1] <- sapply(strsplit(as.character(dat$stand)[idx1], "_"), `[`, 2)
##' idx2 <- seq(2, nrow(dat), by=2)
##' dat$DGE[idx2] <- sapply(strsplit(as.character(dat$stand)[idx2], "_"), `[`, 2)
##' dat$IGE[idx2] <- sapply(strsplit(as.character(dat$stand)[idx2], "_"), `[`, 1)
##' dat$DGE <- factor(dat$DGE, levels=genos)
##' dat$IGE <- factor(dat$IGE, levels=genos)
##' listContr <- list(block="contr.sum")
##' X <- model.matrix(~ 1 + block, data=dat, contrasts=listContr)
##' tmp <- mkZintraspe(dat, colIDfocal="DGE", colIDneighbors="IGE"); Z_D <- tmp$Z_DBV; Z_I <- tmp$Z_SBV
##' truth <- list("intercept"=100,
##'               "var_DGE"=10,
##'               "var_IGE"=2,
##'               "cor_DGE-IGE"=-0.8,
##'               "var_error"=1)
##' set.seed(1234)
##' truth[["blockEffs"]] <- sample(x=c(-1,1), size=nbBlocks - 1, replace=TRUE) *
##'   rnorm(n=nbBlocks - 1, mean=3, sd=5)
##' Sigma <- matrix(c(truth$var_DGE, NA, NA, truth$var_IGE), 2, 2)
##' Sigma[1,2] <- Sigma[2,1] <- truth[["cor_DGE-IGE"]] * sqrt(Sigma[1,1] * Sigma[2,2])
##' tmp <- mvtnorm::rmvnorm(n=nbGenos, sigma=Sigma)
##' truth[["DGEs"]] <- tmp[,1]
##' truth[["IGEs"]] <- tmp[,2]
##' plot(truth[["DGEs"]], truth[["IGEs"]], las=1, main="Simulated genetic effects")
##' abline(h=0, v=0, lty=2)
##' abline(lm(truth[["IGEs"]] ~ truth[["DGEs"]]), col="red")
##' truth[["errors"]] <- rnorm(n=nrow(dat), mean=0, sd=sqrt(truth$var_error))
##' y <- X %*% c(truth$intercept, truth$blockEffs) +
##'   Z_D %*% truth[["DGEs"]] +
##'   Z_I %*% truth[["IGEs"]] +
##'   truth$errors
##' dat$pheno <- y[,1]
##' head(dat)
##' if(FALSE){
##'   hist(dat$pheno, las=1, main="Simulated data")
##'   boxplot(pheno ~ block, data=dat, las=1, main="Simulated data")
##' }
##'
##' ## fit the model
##' myformula <- pheno ~ -1 + block
##' listZ <- list("DGE"=Z_D, "IGE"=Z_I)
##' system.time(
##'    fit <- inlaDGEIGE(myformula, dat, listZ, listContr))
##'
##' ## check the fixed effects' estimates
##' cbind(c(truth$intercept, truth$blockEffs),
##'       fit$summary.random$ID_fix$mean + c(mean(dat$pheno), 0, 0))
##'
##' ## check the (co)variances' estimates
##' fit$summary.hyperpar
##' names(fit$internal.marginals.hyperpar)
##' library(INLA)
##' postSamples <- inla.hyperpar.sample(1000, fit)
##' op <- par(mfrow=c(2,2))
##' hist(postSamples[,1], breaks="FD", main="var_error")
##' abline(v=truth$var_error, col="red")
##' hist(1/postSamples[,2], breaks="FD", main="var_DGE")
##' abline(v=c(0, truth$var_DGE), col="red")
##' hist(1/postSamples[,3], breaks="FD", main="var_IGE")
##' abline(v=c(0, truth$var_IGE), col="red")
##' hist(postSamples[,4], breaks="FD", main="cor_DGE-IGE")
##' abline(v=c(0, truth$`cor_DGE-IGE`), col="red")
##' par(op)
##'
##' ## check the random effects' estimates
##' op <- par(mfrow=c(1,2))
##' idx <- 1:25
##' cor(truth$DGEs, fit$summary.random$ID_iid2d$mean[idx])
##' plot(fit$summary.random$ID_iid2d$mean[idx], truth$DGEs, xlab="BLUPs(DGE)", ylab="true DGEs",
##'      main="Accuracy with INLA", las=1, pch=19)
##' abline(a=0, b=1, v=0, h=0, lty=2)
##' abline(lm(truth$DGEs ~ fit$summary.random$ID_iid2d$mean[idx]), col="red")
##' idx <- (25+1):50
##' cor(truth$IGEs, fit$summary.random$ID_iid2d$mean[idx])
##' plot(fit$summary.random$ID_iid2d$mean[idx], truth$IGEs, xlab="BLUPs(IGE)", ylab="true IGEs",
##'      main="Accuracy with INLA", las=1, pch=19)
##' abline(a=0, b=1, v=0, h=0, lty=2)
##' abline(lm(truth$IGEs ~ fit$summary.random$ID_iid2d$mean[idx]), col="red")
##' par(op)
##' }
##' @noRd
inlaDGEIGE <- function(formula, data, listZ, contrasts = NULL, verbose = FALSE) {
  stopifnot(
    requireNamespace("INLA", quietly = TRUE),
    is.list(listZ),
    all(c("DGE", "IGE") %in% names(listZ))
  )

  ## prepare the matrix of linear predictors
  formFix <- as.character(formula)[3]
  fixEffs <- trimws(strsplit(formFix, "\\+")[[1]])
  fixEffs <- fixEffs[!fixEffs %in% c("1", "-1", "0")]
  if (length(fixEffs) > 0) {
    X <- model.matrix(as.formula(paste0("~ 1 + ", fixEffs)),
      data = data, contrasts = contrasts
    )
    A <- cbind(X, listZ$DGE, listZ$IGE)
    nbFix <- ncol(X)
    nbRndGE <- ncol(listZ$DGE)
    stopifnot(nbFix + 2 * nbRndGE == ncol(A))
    ID_fix <- c(1:nbFix, rep(NA, 2 * nbRndGE))
    ID_iid2d <- c(rep(NA, nbFix), 1:(2 * nbRndGE))
  } else {
    A <- cbind(listZ$DGE, listZ$IGE)
    nbRndGE <- ncol(listZ$DGE)
    stopifnot(2 * nbRndGE == ncol(A))
    ID_iid2d <- 1:(2 * nbRndGE)
  }

  ## prepare the formula
  response <- as.character(formula)[2]
  form <- paste0(response, " ~ -1")
  if (length(fixEffs) > 0) {
    form <- paste0(
      form,
      " + f(ID_fix, model='iid',
                       hyper=list(prec=list(initial=0, fixed=TRUE)))"
    )
  }
  form <- paste0(
    form,
    " + f(ID_iid2d, model='iid2d', n=2*nbRndGE,
                       hyper=list(theta1=list(param=c(4, 0.1, 0.1, 0))))"
  )
  form <- as.formula(form)
  if (verbose) {
    print(form)
  }

  ## prepare the data
  list_dat <- list()
  list_dat[[response]] <- data[[response]]

  ## perform the inference
  fit <- INLA::inla(
    formula = form,
    control.predictor = list(A = A, compute = TRUE, precision = 1e9),
    data = list_dat,
    offset = mean(data[[response]]),
    control.compute = list(dic = TRUE)
  )

  return(fit)
}

##' Fit a GMA-SMA model with BLUPF90
##'
##' Fits a GMA-SMA model with \href{http://nce.ads.uga.edu/wiki/doku.php}{BLUPF90} via the \href{http://famuvie.github.io/breedR/}{breedR} package.
##' @param formFix see \code{\link[breedR:remlf90]{breedR::remlf90()}} with neither \code{"(1|GMA)"} nor \code{"(1|SMA)"}
##' @param data see \code{\link[breedR:remlf90]{breedR::remlf90()}}
##' @param listZ see argument \code{"ZList"} of \code{\link[MM4LMM:MMEst]{MM4LMM::MMEst()}}; the names of this list should contain at least \code{"GMA"} and, optionally, \code{"SMA"}
##' @param listVCov see argument \code{"VarList"} of \code{\link[MM4LMM:MMEst]{MM4LMM::MMEst()}}; the names of this list should be the same as the names of \code{"listZ"}
##' @return list
##' @seealso \code{\link{mkZGMA}}, \code{\link{mkZSMA}}, \code{\link{fitGMASMA}}
##' @author Timothee Flutre
##' @examples
##' \dontrun{
##' ## generate fake data
##' nbGenos <- 25
##' genos <- sprintf("g%02i", 1:nbGenos)
##' pairs <- t(combn(x=genos, m=2))
##' stands <- paste(pairs[,1], pairs[,2], sep="_")
##' nbBlocks <- 3
##' blocks <- LETTERS[1:nbBlocks]
##' dat <- do.call(rbind, lapply(blocks, function(block){
##'   cbind(stands=as.data.frame(stands, stringsAsFactors=TRUE),
##'         block=as.factor(block))
##' }))
##' listContr <- list(block="contr.sum")
##' X <- model.matrix(~ 1 + block, data=dat, contrasts=listContr)
##' Z_GMA <- mkZGMA(df=dat, col="stands", sep="_")
##' truth <- list("intercept"=100,
##'               "var_GMA"=10,
##'               "var_error"=1)
##' set.seed(1234)
##' truth[["blockEffs"]] <- sample(x=c(-1,1), size=nbBlocks - 1, replace=TRUE) *
##'   rnorm(n=nbBlocks - 1, mean=3, sd=5)
##' truth[["GMAs"]] <- rnorm(n=nbGenos, mean=0, sd=sqrt(truth$var_GMA))
##' truth[["errors"]] <- rnorm(n=nrow(dat), mean=0, sd=sqrt(truth$var_error))
##' y <- X %*% c(truth$intercept, truth$blockEffs) +
##'   Z_GMA %*% truth$GMAs +
##'   truth$errors
##' dat$pheno <- y[,1]
##' if(FALSE){
##'   hist(dat$pheno, las=1, main="Simulated data")
##'   boxplot(pheno ~ block, data=dat, las=1, main="Simulated data")
##' }
##'
##' ## fit the model
##' myformFix <- pheno ~ 1 + block
##' listZ <- list("GMA"=Z_GMA)
##' listVCov <- list("GMA"=diag(ncol(Z_GMA)))
##' fit <- bfGMASMA(myformFix, dat, listZ, listVCov)
##'
##' ## check the results
##' BLUEs <- fixef(fit)
##' BLUEs_contr <- solve(cbind(intercept=1, contr.sum(levels(dat$block)))) %*% BLUEs$block
##' data.frame("true"=c(truth$intercept, truth$blockEffs),
##'            "estim"=BLUEs_contr)
##' estV <- fit$var
##' data.frame("true"=c(truth$var_GMA, truth$var_error),
##'            "estim"=fit$var[,1],
##'            row.names=c("GMA","error"))
##' BLUPs <- ranef(fit)
##' cor(truth$GMAs, BLUPs$GMA)
##' plot(BLUPs$GMA, truth$GMAs,
##'      xlab="BLUPs(GMA)", ylab="true GMAs",
##'      main="Accuracy with BLUPF90", las=1, pch=19)
##' abline(a=0, b=1, v=0, h=0, lty=2)
##' abline(lm(truth$GMAs ~ BLUPs$GMA), col="red")
##' }
##' @noRd
bfGMASMA <- function(formFix, data, listZ, listVCov) {
  stopifnot(
    requireNamespace("breedR", quietly = TRUE),
    is.list(listZ),
    !is.null(names(listZ)),
    "GMA" %in% names(listZ),
    is.list(listVCov),
    !is.null(names(listVCov)),
    all(names(listVCov) == names(listZ))
  )

  out <- NULL

  listGeneric <- list(GMA = list(
    listZ[["GMA"]],
    listVCov[["GMA"]]
  ))
  isSMA <- grepl("SMA", names(listZ))
  if (any(isSMA)) {
    for (SMAtype in names(listZ)[which(isSMA)]) {
      listGeneric[[SMAtype]] <- list(
        listZ[[SMAtype]],
        listVCov[[SMAtype]]
      )
    }
  }

  fit <- breedR::remlf90(
    fixed = formFix,
    generic = listGeneric,
    data = data
  )

  return(fit)
}

##' @noRd
tmbGMASMA_prepIn <- function(data, formFix, listZ, listVCov, REML, contrasts) {
  out <- list()

  ## design matrix of fixed effects:
  formChar <- as.character(formFix)
  response <- formChar[2]
  RHS <- formChar[3] # right-hand side
  X <- model.matrix(as.formula(paste("~", RHS)),
    data = data, contrasts.arg = contrasts
  )

  ## design and vcov matrices of random effects:
  listData <- list(
    model = "GMA_SMA",
    y = data[[response]],
    X = X,
    Z_GMA = listZ[["GMA"]],
    K = listVCov[["GMA"]]
  )
  listParams <- list(
    beta = rep(0, ncol(X)),
    GMA = rep(0, ncol(listZ[["GMA"]])),
    log_sigma_GMA = log(1),
    log_sigma_e = log(1)
  )
  vecRnd <- c("GMA")
  if (REML) {
    vecRnd <- c(vecRnd, "beta")
  }

  inMap <- list()
  if (!"SMA" %in% names(listZ) &
    !"SMA_ij" %in% names(listZ) &
    !"SMA_ii" %in% names(listZ)) {
    listData$Z_SMA1 <- matrix(0)
    listData$Z_SMA2 <- matrix(0)
    listParams$SMA1 <- 0
    listParams$log_sigma_SMA1 <- log(1)
    listParams$SMA2 <- 0
    listParams$log_sigma_SMA2 <- log(1)
    inMap[["log_sigma_SMA1"]] <- NA
    inMap[["log_sigma_SMA2"]] <- NA
  }
  if ("SMA" %in% names(listZ)) {
    listData$Z_SMA1 <- listZ[["SMA"]]
    listData$Z_SMA2 <- matrix(0)
    listParams$SMA1 <- rep(0, ncol(listZ[["SMA"]]))
    listParams$log_sigma_SMA1 <- log(1)
    listParams$SMA2 <- 0
    listParams$log_sigma_SMA2 <- log(1)
    vecRnd <- c(vecRnd, "SMA1")
    inMap[["log_sigma_SMA2"]] <- NA
  } else if ("SMA_ij" %in% names(listZ)) {
    listData$Z_SMA1 <- listZ[["SMA_ij"]]
    if ("SMA_ii" %in% names(listZ)) {
      listData$Z_SMA2 <- listZ[["SMA_ii"]]
    }
    listParams$SMA1 <- rep(0, ncol(listZ[["SMA_ij"]]))
    if ("SMA_ii" %in% names(listZ)) {
      listParams$SMA2 <- rep(0, ncol(listZ[["SMA_ii"]]))
    }
    listParams$log_sigma_SMA1 <- log(1)
    if ("SMA_ii" %in% names(listZ)) {
      listParams$log_sigma_SMA2 <- log(1)
    } else {
      inMap[["log_sigma_SMA2"]] <- NA
    }
    vecRnd <- c(vecRnd, "SMA1")
    if ("SMA_ii" %in% names(listZ)) {
      vecRnd <- c(vecRnd, "SMA2")
    }
  }
  inMap <- lapply(inMap, as.factor)

  out$listData <- listData
  out$listParams <- listParams
  out$vecRnd <- vecRnd
  out$map <- inMap

  return(out)
}

##' @noRd
tmbGMASMA_prepOut <- function(outTmb, listData) {
  if (isTRUE(is.na(outTmb$map[["log_sigma_SMA1"]]))) {
    outTmb$report$var_SMA1 <- NULL
    outTmb$report$SMA1 <- NULL
  }
  if (isTRUE(is.na(outTmb$map[["log_sigma_SMA2"]]))) {
    outTmb$report$var_SMA2 <- NULL
    outTmb$report$SMA2 <- NULL
  }
  return(outTmb)
}

##' Fit a GMA-SMA model with TMB
##'
##' Fits a GMA-SMA model with \href{https://cran.r-project.org/package=TMB}{TMB}.
##' @param formFix a symbolic description of the model to be fitted, but without random effects, as an object of class \code{"formula"}
##' @param data data frame containing the variables in the model, but only the fixed effects; see \code{\link[stats:lm]{stats::lm()}}
##' @param listZ named list of design matrices; its names should be \code{"GMA"} (compulsory) and, optionally, one of \code{"SMA"}, \code{"SMA_ij"} or \code{"SMA_ii"}
##' @param listVCov named list of variance-covariance matrices; the names of this list should be the same as the names of \code{"listZ"}; for the moment, the var-cov matrix for the SMAs is not used (TODO: fix this)
##' @param REML logical
##' @param contrasts see \code{\link[stats:model.matrix]{stats::model.matrix()}}
##' @param lOptions named list of options (for experts)
##' @param verbose verbosity level
##' @return list; the "par" component contains the parameters' estimates, "SMA1" corresponding to either "SMA" or "SMA"_ij, and "SMA2" corresponding to "SMA_ii"
##' @seealso \code{\link{mkZGMA}}, \code{\link{mkZSMA}}, \code{\link{lmerGMASMA}}, \code{\link{mmGMASMA}}
##' @author Timothee Flutre
##' @examples
##' \dontrun{
##' ## generate fake data
##' nbGenos <- 25
##' genos <- sprintf("g%02i", 1:nbGenos)
##' pairs <- t(combn(x=genos, m=2))
##' stands <- paste(pairs[,1], pairs[,2], sep="_")
##' nbBlocks <- 3
##' blocks <- LETTERS[1:nbBlocks]
##' dat <- do.call(rbind, lapply(blocks, function(block){
##'   cbind(stands=as.data.frame(stands, stringsAsFactors=TRUE),
##'         block=as.factor(block))
##' }))
##' listContr <- list(block="contr.sum")
##' X <- model.matrix(~ 1 + block, data=dat, contrasts=listContr)
##' Z_GMA <- mkZGMA(df=dat, col="stands", sep="_")
##' truth <- list("intercept"=100,
##'               "var_GMA"=10,
##'               "var_error"=1)
##' set.seed(1234)
##' truth[["blockEffs"]] <- sample(x=c(-1,1), size=nbBlocks - 1, replace=TRUE) *
##'   rnorm(n=nbBlocks - 1, mean=3, sd=5)
##' truth[["GMAs"]] <- rnorm(n=nbGenos, mean=0, sd=sqrt(truth$var_GMA))
##' truth[["errors"]] <- rnorm(n=nrow(dat), mean=0, sd=sqrt(truth$var_error))
##' y <- X %*% c(truth$intercept, truth$blockEffs) +
##'   Z_GMA %*% truth$GMAs +
##'   truth$errors
##' dat$pheno <- y[,1]
##' if(FALSE){
##'   hist(dat$pheno, las=1, main="Simulated data")
##'   boxplot(pheno ~ block, data=dat, las=1, main="Simulated data")
##' }
##'
##' ## fit the model
##' myformFix <- pheno ~ 1 + block
##' listZ <- list("GMA"=Z_GMA)
##' listVCov <- list("GMA"=diag(ncol(Z_GMA)))
##' system.time(
##'   fit <- tmbGMASMA(myformFix, dat, listZ, listVCov, contrasts=listContr))
##'
##' ## check the results
##' BLUEs <- fit$fit$par[1:ncol(X)]
##' data.frame("true"=c(truth$intercept, truth$blockEffs),
##'            "estim"=BLUEs)
##' data.frame("true"=c(truth$var_GMA, truth$var_error),
##'            "estim"=exp(fit$fit$par[paste0("log_sigma_", c("GMA","e"))])^2,
##'            row.names=c("GMA","error"))
##' BLUPs_GMA <- summary(fit$sdrep, select="random")[,1]
##' cor(truth$GMAs, BLUPs_GMA)
##' plot(BLUPs_GMA, truth$GMAs,
##'      xlab="BLUPs(GMA)", ylab="true GMAs",
##'      main="Accuracy with TMB", las=1, pch=19)
##' abline(a=0, b=1, v=0, h=0, lty=2)
##' abline(lm(truth$GMAs ~ BLUPs_GMA), col="red")
##' }
##' @noRd
tmbGMASMA <- function(formFix, data, listZ, listVCov, REML = TRUE,
                      contrasts = NULL, lOptions = NULL, verbose = FALSE) {
  SMAtypes <- c("SMA", "SMA_ij", "SMA_ii")
  withSMAs <- setNames(rep(FALSE, 3), SMAtypes)

  ## quick reformat and check
  stopifnot(
    is.list(listZ),
    !is.null(names(listZ)),
    "GMA" %in% names(listZ),
    is.list(listVCov),
    !is.null(names(listVCov)),
    all(names(listVCov) == names(listZ))
  )
  for (SMAtype in SMAtypes) {
    withSMAs[SMAtype] <- SMAtype %in% names(listZ)
    if (withSMAs[SMAtype]) {
      stopifnot(
        is.matrix(listZ[[SMAtype]]),
        is.matrix(listVCov[[SMAtype]])
      )
    }
  }
  if (!is.null(lOptions)) {
    stopifnot(
      is.list(lOptions),
      !is.null(names(lOptions)),
      all(names(lOptions) %in% c(
        "makeAD", "optim",
        "iter.max", "quantifUncertain"
      ))
    )
  } else {
    lOptions <- list()
  }
  for (step in c("makeAD", "optim", "quantifUncertain")) {
    if (!step %in% names(lOptions)) {
      lOptions[[step]] <- TRUE
    }
  }
  if (!"iter.max" %in% names(lOptions)) {
    lOptions$iter.max <- 200
  }

  out <- list()
  out$REML <- REML

  if (verbose) {
    print("input preparation")
  }
  inputs4TMB <- tmbGMASMA_prepIn(data, formFix, listZ, listVCov, REML, contrasts)
  if (verbose) {
    str(inputs4TMB)
  }

  if (lOptions$makeAD) {
    if (verbose) {
      print("automatic differentiation")
    }
    f <- TMB::MakeADFun(
      data = inputs4TMB$listData,
      parameters = inputs4TMB$listParams,
      map = inputs4TMB$map,
      random = inputs4TMB$vecRnd,
      DLL = "plantmix_TMBExports",
      silent = !verbose
    )
    if (verbose) {
      f$env$tracepar <- TRUE
    } # to debug (it prints param values at each iter)
    out$obj <- f # as in glmmTMB
    out$fn <- f$fn
    out$gr <- f$gr
    out$map <- inputs4TMB$map
  }

  if (lOptions$makeAD & lOptions$optim) {
    if (verbose) {
      print("optimization")
    }
    capture <- capture.output(
      fit <- nlminb(
        start = f$par, objective = f$fn, gradient = f$gr, hessian = NULL,
        control = list("trace" = 1)
      )
    )
    out$fit <- fit
    out$trace <- traceFromNlminb(capture)
    out$report <- f$report() # return derived values; as in glmmTMB
  }

  if (lOptions$makeAD & lOptions$optim & lOptions$quantifUncertain) {
    if (verbose) {
      print("uncertainty quantification")
    }
    out$sdr <- sdreport(f)
    out$sry_sdr <- summary(out$sdr)

    if (verbose) {
      print("output preparation")
    }
    out <- tmbGMASMA_prepOut(out)
  }

  return(out)
}

##' Mixing abilities
##'
##' Extracts general and specific mixing abilities from a merMod object (lme4).
##' @param fit "merMod" object returned by \code{\link[lme4]{lmer}}, as is done by \code{\link{fitGMASMAs}} with \code{pkg="lme4"}
##' @param sep character separating genotype names
##' @return data.frame
##' @seealso \code{\link{fitGMASMA}}
##' @author Timothee Flutre
##' @noRd
getMixingAbilities <- function(fit, sep = "_") {
  out <- NULL

  if (!is(fit, "merMod")) {
    msg <- "does not work (yet) for objects other than merMod (lme4)"
    warning(msg, immediate. = TRUE)
  } else {
    rndEffs <- ranef(fit)
    stopifnot("GMA" %in% names(rndEffs))
    GMAs <- rndEffs$GMA
    withSMA <- "SMA" %in% names(rndEffs)
    withSMA_ii <- "SMA_ii" %in% names(rndEffs)
    withSMA_ij <- "SMA_ij" %in% names(rndEffs)

    if (!withSMA & !withSMA_ii & !withSMA_ij) {
      out <- data.frame(
        geno1 = rownames(GMAs),
        GMA = GMAs[, 1]
      )
    } else if (withSMA & !withSMA_ii & !withSMA_ij) {
      SMAs <- ranef(fit)$SMA
      out <- data.frame(
        geno1 = sapply(strsplit(rownames(SMAs), paste0("\\", sep)), `[`, 1),
        GMA = NA,
        geno2 = sapply(strsplit(rownames(SMAs), paste0("\\", sep)), `[`, 2),
        SMAtype = NA,
        SMA = SMAs[, 1]
      )
      out$GMA <- GMAs[out$geno1, 1]
      out$SMAtype[which(out$geno1 == out$geno2)] <- "SMA_ii"
      out$SMAtype[which(out$geno1 != out$geno2)] <- "SMA_ij"
      idx <- sapply(out, is.character)
      out[idx] <- lapply(out[idx], as.factor)
    } else if (!withSMA & withSMA_ii & withSMA_ij) {
      SMAs_ii <- ranef(fit)$SMA_ii
      stopifnot(all(nrow(SMAs_ii) == nrow(GMAs)))
      SMAs_ij <- ranef(fit)$SMA_ij
      out <- data.frame(
        geno1 = c(
          rownames(GMAs),
          sapply(strsplit(rownames(SMAs_ij), ","), `[`, 1)
        ),
        GMA = NA,
        geno2 = c(
          rownames(GMAs),
          sapply(strsplit(rownames(SMAs_ij), ","), `[`, 2)
        ),
        SMAtype = c(
          rep("SMA_ii", nrow(GMAs)),
          rep("SMA_ij", nrow(SMAs_ij))
        ),
        SMA = c(SMAs_ii[, 1], SMAs_ij[, 1])
      )
      out$GMA <- GMAs[out$geno1, 1]
    } else {
      msg <- "could not detect the model type"
      warning(msg, immediate. = TRUE)
    }

    idx <- sapply(out, is.character)
    out[idx] <- lapply(out[idx], as.factor)
  }

  return(out)
}

## ============================================================================
## code from Emma Forst

## misc.

## Utils for GMA/SMA models
##
## Returns the names of the genotypes in the mixtures of interest.
## @param df data frame
## @param col index or name of the column of \code{df} containing the genotype names per pure stand or mixture
## @param sep character seperating genotype names
## @return character vector
## @author Emma Forst [aut], Timothée Flutre [ctb]
## @examples
## \dontrun{
## dat <- data.frame(stand=rep(paste0("stand", 1:5), 2),
##                   comps=rep(c("geno1-geno1", "geno2-geno2", "geno3-geno3",
##                               "geno1-geno2", "geno2-geno3"),
##                             2),
##                   rep=rep(c("A","B"), each=5),
##                   yield=rnorm(10, 100, 10))
## getMixedGenos(dat, "comps", "-")
## }
getMixedGenos <- function(df, col, sep = ",") {
  sort(unique(unlist(strsplit(as.character(df[, col]), sep))))
}

## Utils for GMA/SMA models
##
## Returns the genotypes of the mixtures of interest.
## @param df data frame
## @param col index or name of the column of \code{df} containing the genotype names per pure stand or mixture
## @param sep character seperating genotype names
## @return list which components are character vectors (one per mixture)
## @author Emma Forst [aut], Timothée Flutre [ctb]
## @examples
## \dontrun{
## dat <- data.frame(stand=rep(paste0("stand", 1:5), 2),
##                   comps=rep(c("geno1-geno1", "geno2-geno2", "geno3-geno3",
##                               "geno1-geno2", "geno2-geno3"),
##                             2),
##                   rep=rep(c("A","B"), each=5),
##                   yield=rnorm(10, 100, 10))
## getGenosPerMix(dat, "comps", "-")
## }
getGenosPerMix <- function(df, col, sep = ",") {
  tmp <- strsplit(as.character(df[, col]), sep)
  unique(lapply(tmp[which(lapply(tmp, length) > 1)], sort))
}

## Utils for GMA/SMA models
##
## Returns pairs of genotypes for all stands (pure and mixed).
## @param mixed.genos output from \code{\link{getMixedGenos}}
## @param genos.per.mix output from \code{\link{getGenosPerMix}}
## @param sep character that will be used to separate the genotype names in a given pair
## @return character vector
## @author Emma Forst [aut]
## @examples
## \dontrun{
## dat <- data.frame(stand=rep(paste0("stand", 1:5), 2),
##                   comps=rep(c("geno1-geno1", "geno2-geno2", "geno3-geno3",
##                               "geno1-geno2", "geno2-geno3"),
##                             2),
##                   rep=rep(c("A","B"), each=5),
##                   yield=rnorm(10, 100, 10))
## mixed_genos <- getMixedGenos(dat, "comps", "-")
## genos_per_mix <- getGenosPerMix(dat, "comps", "-")
## getGenoPairs(mixed_genos, genos_per_mix)
## }
getGenoPairs <- function(mixed.genos, genos.per.mix, sep = "-") {
  Pairs <- lapply(genos.per.mix, function(ll) {
    tmp <- utils::combn(x = ll, m = 2, simplify = FALSE)
    tmp2 <- unlist(lapply(1:length(tmp), function(nn) {
      paste(tmp[[nn]][1], tmp[[nn]][2], sep = sep)
    }))
  })
  pairs.PS <- paste(mixed.genos, mixed.genos, sep = sep) # pure stands
  PairsMixedGenotypes <- sort(c(unique(unlist(Pairs)), pairs.PS))
  return(PairsMixedGenotypes)
}

## Utils for GMA/SMA models
##
## Utils
## @param df data frame
## @param response.names character vector with the names of responses
## @return list which components are vectors containing phenotypic values (one per trait)
## @author Emma Forst [aut], Timothée Flutre [ctb]
## @examples
## \dontrun{
## dat <- data.frame(stand=rep(paste0("stand", 1:5), 2),
##                   comps=rep(c("geno1-geno1", "geno2-geno2", "geno3-geno3",
##                               "geno1-geno2", "geno2-geno3"),
##                             2),
##                   rep=rep(c("A","B"), each=5),
##                   yield=rnorm(10, 100, 10))
## getPhenos(dat, "yield")
## }
getPhenos <- function(df, response.names) {
  out <- lapply(response.names, function(pp) {
    Reduce("c", df[[pp]])
  })
  names(out) <- response.names
  return(out)
}

## Design matrix of fixed effects for GMA/SMA models
##
## Design matrix for fixed effects.
## @param df data frame
## @param ListX character vector with the names of fixed effects
## @return matrix
## @author Emma Forst [aut], Timothée Flutre [ctb]
## @examples
## \dontrun{
## dat <- data.frame(stand=rep(paste0("stand", 1:5), 2),
##                   comps=rep(c("geno1-geno1", "geno2-geno2", "geno3-geno3",
##                               "geno1-geno2", "geno2-geno3"),
##                             2),
##                   rep=rep(c("A","B"), each=5),
##                   yield=rnorm(10, 100, 10))
## (X <- BuildX(dat))
## (X <- BuildX(dat, "rep"))
## }
BuildX <- function(df, ListX = NULL) {
  if (is.null(ListX)) {
    X <- as.matrix(rep(1, nrow(df)))
  } else {
    Var <- Reduce(
      "rbind",
      df[, colnames(df) %in% ListX, drop = F]
    )
    Fake <- data.frame(Y = stats::rnorm(length(Var)), Var)
    X <- stats::model.matrix(stats::lm(
      formula = Y ~ .,
      data = Fake,
      contrasts = list(Var = "contr.sum")
    ))
    colnames(X)[-1] <- ListX
  }
  return(X)
}

## Design matrix of GMAs for GMA/SMA models
##
## Design matrix for GMAs.
## @param df data frame
## @param Gcol index of the column containing the genotypes per pure stand or mixture
## @param G character vector of genotype names (output from \code{\link{getMixedGenos}})
## @param sep separator
## @return matrix
## @author Emma Forst [aut]
## @examples
## \dontrun{
## dat <- data.frame(stand=rep(paste0("stand", 1:5), 2),
##                   comps=rep(c("geno1-geno1", "geno2-geno2", "geno3-geno3",
##                               "geno1-geno2", "geno2-geno3"),
##                             2),
##                   rep=rep(c("A","B"), each=5),
##                   yield=rnorm(10, 100, 10))
## mixed_genos <- getMixedGenos(dat, "comps", "-")
## (Z_GMA <- BuildZ_GMA(dat, 2, mixed_genos, "-"))
## }
BuildZ_GMA <- function(df, Gcol, G, sep = ",") {
  Z <- Reduce("rbind", lapply(1:nrow(df), function(rr) {
    MixtSize <- length(strsplit(as.character(df[rr, Gcol]), sep)[[1]])
    Reduce("+", lapply(1:MixtSize, function(kk) {
      genolignetmp <- strsplit(as.character(df[, Gcol]), sep)[[rr]][[kk]]
      sapply(G, function(gg) {
        genolignetmp == gg
      })
    })) / MixtSize
  }))
  rownames(Z) <- NULL
  return(Z)
}

## Design matrix of SMAs for GMA/SMA models
##
## Design matrix for SMAs (model 2).
## @param df data frame
## @param Gcol index of the column containing the genotypes per pure stand or mixture
## @param G character vector of genotype names (output from \code{\link{getMixedGenos}})
## @param Gpairs character vector of genotype pairs (output from \code{\link{getGenoPairs}})
## @param sep separator
## @return matrix
## @author Emma Forst [aut]
## @examples
## \dontrun{
## dat <- data.frame(stand=rep(paste0("stand", 1:5), 2),
##                   comps=rep(c("geno1-geno1", "geno2-geno2", "geno3-geno3",
##                               "geno1-geno2", "geno2-geno3"),
##                             2),
##                   rep=rep(c("A","B"), each=5),
##                   yield=rnorm(10, 100, 10))
## mixed_genos <- getMixedGenos(dat, "comps", "-")
## genos_per_mix <- getGenosPerMix(dat, "comps", "-")
## geno_pairs <- getGenoPairs(mixed_genos, genos_per_mix)
## (Z_SMA <- BuildZ_SMAmodel2(dat, 2, mixed_genos, geno_pairs, "-"))
## }
BuildZ_SMAmodel2 <- function(df, Gcol, G, Gpairs, sep = ",") {
  Z <- Reduce("rbind", lapply(1:nrow(df), function(rr) {
    MixtSize <- length(strsplit(as.character(df[rr, Gcol]), sep)[[1]])
    if (MixtSize > 1) {
      genolignetmp <- strsplit(as.character(df[, Gcol]), sep)[[rr]]
      Comb <- utils::combn(x = genolignetmp, simplify = T, m = 2)
      tmp <- Reduce("+", lapply(1:ncol(Comb), function(CC) {
        sapply(Gpairs, function(gg) {
          paste(pmin(Comb[1, CC], Comb[2, CC]), pmax(Comb[1, CC], Comb[2, CC]), sep = "-") == gg
        })
      })) / choose(MixtSize, 2)
    } else {
      ## no mixture
      tmp <- sapply(Gpairs, function(gg) {
        paste(df[rr, Gcol], df[rr, Gcol], sep = "-") == gg
      })
    }
  }))

  rownames(Z) <- NULL
  return(Z)
}

## Design matrix of SMAs for GMA/SMA models
##
## Design matrix for SMAs (model 3).
## @param df data frame
## @param Gcol index of the column containing the genotypes per pure stand or mixture
## @param Gpairs character vector of genotype pairs (output from \code{\link{getGenoPairs}})
## @param sep separator
## @return matrix
## @author Emma Forst [aut]
## @examples
## \dontrun{
## dat <- data.frame(stand=rep(paste0("stand", 1:5), 2),
##                   comps=rep(c("geno1-geno1", "geno2-geno2", "geno3-geno3",
##                               "geno1-geno2", "geno2-geno3"),
##                             2),
##                   rep=rep(c("A","B"), each=5),
##                   yield=rnorm(10, 100, 10))
## mixed_genos <- getMixedGenos(dat, "comps", "-")
## genos_per_mix <- getGenosPerMix(dat, "comps", "-")
## geno_pairs <- getGenoPairs(mixed_genos, genos_per_mix)
## (Z_SMA <- BuildZ_SMAmodel3(dat, 2, geno_pairs, "-"))
## }
BuildZ_SMAmodel3 <- function(df, Gcol, Gpairs, sep = ",") {
  Z <- Reduce("rbind", lapply(1:nrow(df), function(rr) {
    MixtSize <- length(strsplit(as.character(df[rr, Gcol]), sep)[[1]])
    if (MixtSize > 1) {
      genolignetmp <- strsplit(as.character(df[, Gcol]), sep)[[rr]]
      Comb <- utils::combn(x = genolignetmp, simplify = T, m = 2)
      tmp <- Reduce("+", lapply(1:ncol(Comb), function(CC) {
        sapply(Gpairs, function(gg) {
          paste(pmin(Comb[1, CC], Comb[2, CC]), pmax(Comb[1, CC], Comb[2, CC]), sep = "-") == gg
        })
      })) * (2 / MixtSize**2)

      ## intra-genotypic interactions (i = j)
      tmp <- tmp + Reduce("+", lapply(1:MixtSize, function(kk) {
        sapply(Gpairs, function(gg) {
          paste(genolignetmp[[kk]], genolignetmp[[kk]], sep = "-") == gg
        })
      })) / (MixtSize**2)
    } else {
      ## no mixture
      tmp <- sapply(Gpairs, function(gg) {
        paste(df[rr, Gcol], df[rr, Gcol], sep = "-") == gg
      })
    }
  }))

  rownames(Z) <- NULL
  return(Z)
}

## ============================================================================

## lme4 hacking

## lme4 hacking for GMA-SMA models
##
## Parse the data and formula + information on random effects structure.
## @param ListZ list of design matrices for random effects
## @return list
## @author Emma Forst [aut]
myReTrms <- function(ListZ) {
  reTrms <- list()

  ## fusionne les colonnes de ListZ, prend sa transposee et cree une matrice creuse
  reTrms$Zt <- Matrix::Matrix(t(Reduce("cbind", ListZ)), sparse = TRUE)

  ## Parametrization (numeric vector of variance component parameters) : Initial Value of the covariance parameters
  reTrms$theta <- rep(1, length(ListZ))

  ## an integer vector of indices (of the same length as the X slot at the lambda slot) determining the mapping of the elements of the theta vector to the "x" slot of Lambdat
  reTrms$Lind <- rep(1:length(ListZ), unlist(lapply(ListZ, ncol)))
  reTrms$Gp <- as.integer(unname(cumsum(c(0, unlist(lapply(ListZ, ncol))))))

  ## lower bounds on the covariance parameters
  reTrms$lower <- rep(0, length(ListZ))

  ## t Covariance factor : transpose of the sparse relative covariance factor (or lower triangular relative variance)
  reTrms$Lambdat <- Matrix::Matrix(diag(rep(1, sum(unlist(lapply(ListZ, ncol))))),
    sparse = TRUE
  )

  ## obtient une liste avec les transposees des matrices indicatrices
  reTrms$Ztlist <- lapply(ListZ, function(Z) Matrix::Matrix(t(Z), sparse = TRUE))

  ## les noms des elements de ListZ
  reTrms$cnms <- as.list(names(ListZ))
  names(reTrms$cnms) <- names(ListZ)

  ## list of grouping factors used in the random-effects terms
  ## Flist is Not very clean (to say the least... )
  reTrms$flist <- lapply(ListZ, function(Z) {
    flist.factor <- as.factor(colnames(Z)[apply(Z, 1, function(x) {
      which(stats::rmultinom(n = 1, size = 1, prob = x) == 1)
    })])
    levels(flist.factor) <- colnames(Z)
    return(flist.factor)
  })

  return(reTrms)
}

## lme4 hacking for GMA/SMA models
##
## Fit the model.
## @param Response numeric vector with phenotypic values
## @param X design matrix for fixed effects; see \code{\link{BuildX}}
## @param ListZ list of design matrices for random effects; see \code{\link{BuildZ_GMA}}, \code{\link{BuildZ_SMAmodel2}} and \code{\link{BuildZ_SMAmodel3}}
## @param REML logical
## @param boundary.tol boundary tolerance
## @return merMod object
## @author Emma Forst [aut]
## @seealso \code{\link{ranefGMASMA_EF}}
lmerGMASMA_EF <- function(Response, X, ListZ, REML = TRUE, boundary.tol = 10e-10) {
  notNA <- !is.na(Response)
  Response <- Response[notNA]
  X <- X[notNA, , drop = FALSE]
  ListZ <- lapply(ListZ, function(Z) Z <- Z[notNA, ])

  fr <- stats::model.frame(
    formula = Response ~ .,
    data.frame(Response = Response, X)
  )

  reTrms <- myReTrms(ListZ)

  ##  Create the deviance function (objective function) to be optimized (minimized) to estimate the parameters : return a function to calculate deviance (or restricted deviance) as a function over theta (the random-effect parameters)
  devfun <- lme4::mkLmerDevfun(fr, X, reTrms,
    REML = REML,
    boundary.tol = boundary.tol
  )

  ## pour optimiser les resultats (return the results of an optmization of the deviance function) : character - name of optimizing function(s). A list of functions : "bobyqa"  (from the minqa package) with two par (best-fit parameters), a fval (best-fit function value) and number of function evaluations: 33 ... see running opt[]
  opt <- lme4::optimizeLmer(devfun, boundary.tol = 10e-10, optimizer = "bobyqa")

  return(lme4::mkMerMod(rho = environment(devfun), opt = opt, reTrms = reTrms, fr = fr))
}

## lme4 hacking for GMA/SMA models
##
## Extract the BLUPs (conditional modes of the random effects).
## @param model output from \code{\link{lmerGMASMA_EF}}
## @param condVar logical
## @return ranef.mer object
## @author Emma Forst [aut]
ranefGMASMA_EF <- function(model, condVar = TRUE) { # conditional variance: condVar : indicates if the conditional variance-covariance matrices of the random effects given the response should be added as an attribute (to get the precision of the estimates).
  # formula of the variances of U|Y : D - D Zt solve(Sigma) Z D with Sigma = Z D Zt + S and S=sigma^2 the residual variance
  # re.cond.mode<-tapply(model@u,model@pp$Lind,function(x) x) # u is the vector of conditional model of spherical random effects coefficients
  # names(re.cond.mode)<- names(model@cnms)
  ## Lind is an integer vector that maps the elements of theta to the non-zero in Lambda
  re.cond.mode <- sapply(1:max(model@pp$Lind), function(i) {
    u <- model@u[model@pp$Lind == i]
    Lambdat <- model@pp$Lambdat[model@pp$Lind == i, model@pp$Lind == i]
    b <- as.numeric(t(Lambdat) %*% u)
  })
  if (is.numeric(re.cond.mode)) re.cond.mode <- list(re.cond.mode)
  names(re.cond.mode) <- names(model@cnms)

  if (condVar) {
    Zt <- model@pp$Zt
    D <- stats::sigma(model)^2 * t(model@pp$Lambdat) %*% model@pp$Lambdat # matrice de la variance des blups (sur la diagonale)
    Sigma <- t(Zt) %*% D %*% Zt + stats::sigma(model)^2 * diag(rep(1, ncol(Zt))) # sigma(model) is the variance of residuals
    var.cond <- D - D %*% Zt %*% solve(Sigma) %*% t(Zt) %*% D
    var.cond <- diag(var.cond)
    var.cond.mode <- tapply(var.cond, model@pp$Lind, function(x) x) # cree une liste pour separer les variances des GMA et des SMA
  }

  for (i in 1:length(re.cond.mode)) {
    re.cond.mode[[i]] <- data.frame(re.cond.mode[i])
    names(re.cond.mode[[i]]) <- names(re.cond.mode[i])
    row.names(re.cond.mode[[i]]) <- levels(model@flist[[i]])

    if (condVar) attr(re.cond.mode[[i]], "postVar") <- array(var.cond.mode[[i]], c(1, 1, nrow(re.cond.mode[[i]])))
  }
  attr(re.cond.mode, "class") <- "ranef.mer"
  re.cond.mode
}

##' Proportion of variance explained for GMA-SMA models
##'
##' Calculating the R^2 of the models for each response variable using the formula of equation 1 in Alexander, Tropsha and Winkler (2015):
##' R^2 = 1 - (sum(squared(yobs - ypred)) / sum(squared(yobs - ymean)))
##' @param fit merMod object
##' @return numeric
##' @author Emma Forst [aut]
##' @noRd
myRsq <- function(fit) {
  1 - (sum((fit@resp$y - fit@resp$mu)^2) /
    sum((fit@resp$y - fit@beta[1])^2))
}

## ============================================================================

## reformat lme4 outputs

## ##' Reformat
## ##'
## ##' Reformat.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_Beta_lme4<-function(PhenoList){
##   Beta<-data.frame()
##   Betatmp<-vector("list",length(PhenoList))
##   Betatmp<-lapply(PhenoList, function(pheno,i){
##     i=match(pheno,PhenoList)
##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)
##     btmp1<-summary(mymod.gma)$coefficients
##     btmp2<-summary(mymod)$coefficients
##     btmp3<-summary(mymod.old)$coefficients
##     colnames(btmp1)[[1]]<-paste("Estimate","mymod.gma", sep="_")
##     colnames(btmp2)[[1]]<-paste("Estimate","mymod", sep="_")
##     colnames(btmp3)[[1]]<-paste("Estimate","mymod.old", sep="_")
##     Betatmp[[i]]<-data.frame(phenotype=rep(PhenoList[[i]], length(FixedVarList)+1), fixed=c("mu", FixedVarList),cbind(btmp1,btmp2,btmp3))
##   })
##   Beta<-Reduce('rbind',Betatmp)
## write.csv2(Beta, file="beta.csv", row.names = FALSE)
## }

## ##' Reformat
## ##'
## ##' Reformat.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_VarComp_lme4<-function(PhenoList){
##   VarComp<-data.frame()
##   VarComptmp<-vector("list",length(PhenoList))
##   VarComptmp<-lapply(PhenoList, function(pheno,i){
##     i=match(pheno,PhenoList)

##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)

##     vctmp1<-as.data.frame(VarCorr(mymod.gma),comp=c("Variance","Std.Dev."))[,c(2,4,5)]
##     vctmp1$var1[2]<-"Residual"
##     vctmp1<- rbind(vctmp1[1,],rep(NA,3),vctmp1[-1,])
##     colnames(vctmp1)<-c("mymod.gma","Variance","Std.Dev.")
##     vctmp2<-as.data.frame(VarCorr(mymod),comp=c("Variance","Std.Dev."))[,c(2,4,5)]
##     vctmp2$var1[3]<-"Residual"
##     colnames(vctmp2)<-c("mymod","Variance","Std.Dev.")
##     vctmp3<-as.data.frame(VarCorr(mymod.old),comp=c("Variance","Std.Dev."))[,c(2,4,5)]
##     vctmp3$var1[3]<-"Residual"
##     colnames(vctmp3)<-c("mymod.old","Variance","Std.Dev.")
##     vctmp<-cbind(vctmp1,vctmp2,vctmp3)

##     VarComptmp[[i]]<-data.frame(phenotype=PhenoList[[i]], vctmp)
##   })
##   VarComp<-Reduce('rbind',VarComptmp)
##   write.csv2(VarComp, file="varComp.csv")
## }

## ##' Reformat
## ##'
## ##' Reformat.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_GMA_lme4<-function(PhenoList){
##   GMAtmp<-vector("list",length(PhenoList))
##   GMAtmp<-lapply(PhenoList, function(pheno,i){
##     i=match(pheno,PhenoList)
##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)
##     GMAtmp[[i]]<-data.frame(myranef(mymod.gma)$GMA,myranef(mymod)$GMA, myranef(mymod.old)$GMA, "conditional variances"=attributes(myranef(mymod.gma)$GMA)$postVar[,,],"conditional variances"=attributes(myranef(mymod)$GMA)$postVar[,,],"conditional variances"=attributes(myranef(mymod.old)$GMA)$postVar[,,])
##     colnames(GMAtmp[[i]])<-c("GMA_mymod.gma","GMA_mymod","GMA_mymod.old",rep("conditional.variances",3))
##     write.csv2(GMAtmp[[i]], file=paste(paste("gma",PhenoList[[i]],sep="_"),".csv",sep=""))
##   })
## }

## ##' Reformat
## ##'
## ##' Reformat.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_SMA_lme4<-function(PhenoList){
##   SMAtmp<-vector("list",length(PhenoList))
##   SMAtmp<-lapply(PhenoList, function(pheno,i){
##     i=match(pheno,PhenoList)
##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)
##     SMAtmp[[i]]<-data.frame(myranef(mymod)$SMA, myranef(mymod.old)$SMA, "conditional variances"=attributes(myranef(mymod)$SMA)$postVar[,,],"conditional variances"=attributes(myranef(mymod.old)$SMA)$postVar[,,])
##     colnames(SMAtmp[[i]])<-c("SMA_mymod","SMA_mymod.old",rep("conditional.variances",2))
##     write.csv2(SMAtmp[[i]], file=paste(paste("sma",PhenoList[[i]],sep="_"),".csv",sep=""))
##   })
## }

## ##' Reformat
## ##'
## ##' Reformat.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_Ypred_lme4<-function(PhenoList){
##   Ypredtmp<-vector("list",length(PhenoList))
##   Ypredtmp<-lapply(PhenoList, function(pheno,i){
##     i=match(pheno,PhenoList)
##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)
##     Ypred.mod <- X%*%mymod@beta + Reduce('+',lapply(1:length(ListZ), function(zz){
##       ListZ[[zz]]%*%data.matrix(myranef(mymod)[[zz]])
##       }))
##     Ypred.gma <- X%*%mymod.gma@beta + ListZ[[1]]%*%data.matrix(myranef(mymod.gma)[[1]])
##     Ypred.old <- X%*%mymod.old@beta + Reduce('+',lapply(1:length(ListOldZ), function(zz){
##       ListOldZ[[zz]]%*%data.matrix(myranef(mymod.old)[[zz]])
##       }))
##     Ypredtmp[[i]]<-cbind(Phenotype[[i]],Ypred.gma, Ypred.mod,Ypred.old)
##     colnames(Ypredtmp[[i]])=c("Yobs","mod.gma", "mod", "mod.old")
##     write.csv2(Ypredtmp[[i]], file=paste(paste("ypred",PhenoList[[i]],sep="_"),".csv",sep=""))
##     })
## }

## ##' Reformat
## ##'
## ##' Save plots of observed and predicted values.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_Ypred_plots_lme4<-function(PhenoList){
##   for (i in 1:length(PhenoList)){
##     mymod<- mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListZ = ListOldZ,REML=TRUE)
##     Ypred <- X%*%mymod@beta + Reduce('+',lapply(1:length(ListZ), function(zz){
##       ListZ[[zz]]%*%data.matrix(myranef(mymod)[[zz]])
##     }))
##     Ypred.gma <- X%*%mymod.gma@beta + ListZ[[1]]%*%data.matrix(myranef(mymod.gma)[[1]])
##     Ypred.old <- X%*%mymod.old@beta + Reduce('+',lapply(1:length(ListOldZ), function(zz){
##       ListOldZ[[zz]]%*%data.matrix(myranef(mymod.old)[[zz]])
##     }))
##     tiff(file=paste("Yobs-Ypred_",PhenoList[[i]],".tiff",sep=""),width=3000,height=2500, res=350)
##     par(mar=c(5,4,2,2)+0.1)
##     plot(Phenotype[[i]],Ypred.gma,pch=16,xlab='Observed values',ylab='Predicted values for observations',main = PhenoList[[i]],
##          xlim=c(min(Phenotype[[i]],Ypred.gma,Ypred,Ypred.old, na.rm=TRUE),max(Phenotype[[i]],Ypred.gma,Ypred,Ypred.old, na.rm=TRUE)),
##          ylim=c(min(Phenotype[[i]],Ypred.gma,Ypred,Ypred.old, na.rm=TRUE),max(Phenotype[[i]],Ypred.gma,Ypred,Ypred.old, na.rm=TRUE)))
##     points(Phenotype[[i]],Ypred,col=2,pch=16)
##     points(Phenotype[[i]],Ypred.old,col=3,pch=16)
##     abline(0,1,col=1,lwd=2)
##     Legend <- paste0('RMSE(Y,',c('GMA','SMA','OM'),')=',round(c(RMSE.merMod(mymod.gma),RMSE.merMod(mymod),RMSE.merMod(mymod.old)),3))
##     legend("bottomright",legend=Legend,col=1:3,pch=16)
##     dev.off()
##   }
## }

## ##' Reformat
## ##'
## ##' Save plots of observed and predicted values : en inversant les axes des x et des y.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_Ypred_plots_lme4_formated<-function(PhenoList){
##   PhenoNames<-c("Yield","Spike density", "Spike number per plant", "Grain number per spike", "TKW")
##   for (i in 1:length(PhenoList)){
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod<- mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListZ = ListOldZ,REML=TRUE)
##     Ypred.gma <- mymod.gma@resp$mu
##     Ypred.old <- mymod.old@resp$mu
##     Ypred <- mymod@resp$mu
##     y<-mymod@resp$y
##     tiff(file=paste("Yobs-Ypred_",PhenoList[[i]],".tiff",sep=""),width=14.2,height=12.7, units="cm", res=350)
##     par(mar=c(4.5,4,2.5,1.5)+0.1)
##     plot(Ypred.gma, y,pch=16,cex=0.8,xlab='Predicted values for observations',ylab='Observed values',main = PhenoNames[[i]],ps=10,
##         xlim=c(min(y,Ypred.gma,Ypred,Ypred.old, na.rm=TRUE),max(y,Ypred.gma,Ypred,Ypred.old, na.rm=TRUE)),
##         ylim=c(min(y,Ypred.gma,Ypred,Ypred.old, na.rm=TRUE),max(y,Ypred.gma,Ypred,Ypred.old, na.rm=TRUE)))
##     points(Ypred.old,y,col=2,pch=16,cex=0.8)
##     points(Ypred,y,col="dodgerblue",pch=16, cex=0.8)
##     abline(0,1,col=1,lwd=1)
##     Legend <- paste0('RMSE(',c('Model 1','Model 2','Model 3'),',Y)=',round(c(RMSE.merMod(mymod.gma),RMSE.merMod(mymod.old),RMSE.merMod(mymod)),3))
##     legend("bottomright",legend=Legend,col=c(1,2,"dodgerblue"),pch=16, cex=0.9, pt.cex = 0.8)
##     dev.off()
##   }
## }

## ##' Reformat
## ##'
## ##' Save confidence intervals for fixed effects and variance components.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_confint_lme4<- function(PhenoList){
##   confinterv<-list()
##   for (i in 1:length(PhenoList)) {
##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)
##     confinterv[[i]]<-list(stats::confint(stats::profile(mymod.gma),level=0.9),
##                           stats::confint(stats::profile(mymod),level=0.9),
##                           stats::confint(stats::profile(mymod.old),level=0.9))
##     names(confinterv[[i]])<-c("mymod.gma","mymod","mymod.old")
##   }
##   names(confinterv)<-PhenoList
##   capture.output(confinterv, file="confidence-intervals.txt")
## }

## ##' Reformat
## ##'
## ##' Save Anova symmary for model comparisons.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_model_comparison_lme4<- function(PhenoList){
##   modcomparison<-list()
##   for (i in 1:length(PhenoList)){
##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)
##     modcomparison[[i]]<-list(data.frame(stats::anova(mymod.gma,mymod)),
##                              data.frame(stats::anova(mymod.gma, mymod.old)))
##   }
##   names(modcomparison)<-PhenoList
##   capture.output(modcomparison, file="Models-comparison.txt")
## }

## ##' Reformat
## ##'
## ##' Save model accuracy values.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_model_accuracy_lme4<-function(PhenoList){
##   modaccuracy<-list()
##   for (i in 1:length(PhenoList)){
##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)
##     modaccuracy[[i]]<-as.data.frame(t(data.frame(summary(mymod.gma)$AICtab, summary(mymod)$AICtab, summary(mymod.old)$AICtab))[,-5])
##     # the df of the residual are removed from the data.frame
##     rownames(modaccuracy[[i]])=c("mymod.gma","mymod","mymod.old")
##     modaccuracy[[i]]$RMSE<-c(RMSE.merMod(mymod.gma),RMSE.merMod(mymod),RMSE.merMod(mymod.old))
##   }
##   names(modaccuracy)<-PhenoList
##   capture.output(modaccuracy, file="Models-accuracy.txt")
## }

## ##' Reformat
## ##'
## ##' Save model accuracy values.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_model_accuracy_lme4REML<-function(PhenoList){
##   modaccuracy<-list()
##   for (i in 1:length(PhenoList)){
##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)
##     modaccuracy[[i]]<-as.data.frame(t(data.frame(summary(mymod.gma)$AICtab, summary(mymod)$AICtab, summary(mymod.old)$AICtab)))
##     # the df of the residual are removed from the data.frame
##     rownames(modaccuracy[[i]])=c("mymod.gma","mymod","mymod.old")
##     modaccuracy[[i]]$RMSE<-c(RMSE.merMod(mymod.gma),RMSE.merMod(mymod),RMSE.merMod(mymod.old))
##   }
##   names(modaccuracy)<-PhenoList
##   capture.output(modaccuracy, file="Models-accuracy.txt")
## }

## ##' Reformat
## ##'
## ##' Save model accuracy values.
## ##' @param PhenoList character vector with trait names
## ##' @return nothing
## ##' @author Emma Forst [aut]
## ##' @noRd
## Save_model_accuracy_lme4REML2<-function(PhenoList){
##   modaccuracy<-list()
##   for (i in 1:length(PhenoList)){
##     mymod <-mylmer(Response =  Phenotype[[i]],X, ListZ,REML=TRUE)
##     mymod.gma<-mylmer(Response =  Phenotype[[i]],X, ListZ[1],REML=TRUE)
##     mymod.old<-mylmer(Response =  Phenotype[[i]],X, ListOldZ,REML=TRUE)
##     modaccuracy[[i]]<-as.data.frame(t(data.frame(summary(mymod.gma)$AICtab, summary(mymod)$AICtab, summary(mymod.old)$AICtab)))
##     # the df of the residual are removed from the data.frame
##     rownames(modaccuracy[[i]])=c("mymod.gma","mymod","mymod.old")
##     modaccuracy[[i]]$RMSE<-c(RMSE.merMod(mymod.gma),RMSE.merMod(mymod),RMSE.merMod(mymod.old))
##     modaccuracy[[i]]$R2<-c(myRsq(model=mymod.gma),myRsq(model=mymod),myRsq(model=mymod.old))
##   }
##   names(modaccuracy)<-PhenoList
##   capture.output(modaccuracy, file="Models-accuracy.txt")
## }
