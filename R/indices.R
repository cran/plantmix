## Functions implementing indices useful to analyze mixture data.

##' Relative yield per plant (RYP)
##'
##' Computes the relative yield per plant (RYP) as defined in Fowler (1982, \doi{10.2307/2259865}): RYP_ij = Y_ij / (p Y_i) where Y_ij is the yield of i when mixed with j at a proportion of p, and Y_i is the yield of i in monovarietal culture.
##' Note that, unfortunately, the equation of RYP was called "relative yield" (RY) in table 2 of Williams and McCarthy (2001, \doi{10.1046/j.1440-1703.2001.00368.x}), although these authors do have the right equation, i.e., they only dropped the "per plant".
##' Missing values are forbidden.
##' @param dat data.frame with one yield value per row and with at least four columns: the identifier of each stand, the identifier of the focal genotype or species in each stand, the sowing (plant) proportion of the focal in each stand, and the focal yield in each stand (see the example below); a stand with a single proportion equals to 1 will be interpreted as a monovarietal culture ("pure stand"), a mixture otherwise (and its proportions should sum to 1); note that contrary to \code{\link{RY}}, the proportions will be used to compute the values of RYP
##' @param colIDstand column name for the stand identifiers
##' @param colIDfocal column name for the focal identifiers
##' @param colProp column name for the proportions
##' @param colY column name for the yield values
##' @return input data.frame with an extra column named "RYP"
##' @author Timothee Flutre
##' @seealso \code{\link{RY}}, \code{\link{RYT}}, \code{\link{RYM}}
##' @examples
##' (dat <- data.frame(ID=c("geno1", "geno2", "mixg1g2", "mixg1g2"),
##'                    focal=c("geno1", "geno2", "geno1", "geno2"),
##'                    prop=c(1, 1, 0.5, 0.5),
##'                    yield=c(50, 40, 25, 22)))
##' RYP(dat)
##' @export
RYP <- function(dat, colIDstand = "ID", colIDfocal = "focal", colProp = "prop",
                colY = "yield") {
  stopifnot(
    is.data.frame(dat),
    all(c(colIDstand, colIDfocal, colProp, colY) %in% colnames(dat)),
    is.numeric(dat[[colProp]]),
    is.numeric(dat[[colY]]),
    all(!is.na(dat[[colProp]])),
    all(!is.na(dat[[colY]]))
  )

  out <- data.frame(dat,
    RYP = NA
  )

  ## monovarietal cultures
  isMono <- sapply(dat[[colProp]], function(x){
    isTRUE(all.equal(1.0, x))
  })
  stopifnot(
    sum(isMono) > 0,
    sum(isMono) < nrow(dat)
  )
  out$RYP[which(isMono)] <- 1.0

  ## mixtures
  idxMixes <- which(!isMono)
  mixSumProps <- tapply(
    dat[[colProp]][idxMixes],
    as.character(dat[[colIDstand]][idxMixes]),
    sum
  )
  tmp <- Map(function(x) {
    all.equal(1.0, x)
  }, mixSumProps)
  tmp <- unlist(lapply(tmp, function(x) {
    if (is.logical(x)) {
      x
    } else {
      FALSE
    }
  }))
  if (!all(tmp)) {
    msg <- "not all props are equal to 1"
    stop(msg)
  }
  for (i in idxMixes) {
    mixYield <- dat[[colY]][i]
    prop <- dat[[colProp]][i]
    focalID <- as.character(dat[[colIDfocal]][i])
    idx <- which(dat[[colIDstand]] == focalID)
    stopifnot(length(idx) == 1)
    monoYield <- dat[[colY]][idx]
    out$RYP[i] <- mixYield / (prop * monoYield)
  }

  return(out)
}

##' Relative yield (RY)
##'
##' Computes the "relative yield" (RY) as defined in Fowler (1982, \doi{10.2307/2259865}) citing \href{https://research.wur.nl/en/publications/on-competition/}{de Wit (1960)}: RY_ij = Y_ij / Y_i where Y_ij is the yield of i when mixed with j and Y_i is the yield of i in monovarietal culture.
##' Note that equation 1 of table 2 in Williams and McCarthy (2001, \doi{10.1046/j.1440-1703.2001.00368.x}) is called "RY" but in fact it corresponds to "RYP" (relative yield per plant).
##' Note also that the equation for "RY" in Reiss and Drinkwater (2018, \doi{10.1002/eap.1629}) in fact corresponds to "RYM" (relative yield of mixture).
##' Missing values are forbidden.
##' @param dat data.frame with one yield value per row and with at least four columns: the identifier of each stand, the identifier of the focal genotype or species in each stand, the sowing (plant) proportion of the focal in each stand, and the focal yield in each stand (see the example below); a stand with a single proportion equals to 1 will be interpreted as a monovarietal culture ("pure stand") and a mixture otherwise; note that contrary to \code{\link{RYP}}, the proportions will not be used for anything else
##' @param colIDstand column name for the stand identifiers
##' @param colIDfocal column name for the focal identifiers
##' @param colProp column name for the proportions
##' @param colY column name for the yield values
##' @param avgIfReps if TRUE, replicates of monovarietals will be averaged; if FALSE, the presence of replicates will return an error
##' @return input data.frame with an extra column named "RY"
##' @author Timothee Flutre
##' @seealso \code{\link{estimRYRep}}, \code{\link{RYP}}, \code{\link{RYT}}, \code{\link{RYM}}
##' @examples
##' (dat <- data.frame(ID=c("geno1", "geno2", "mixg1g2", "mixg1g2"),
##'                    focal=c("geno1", "geno2", "geno1", "geno2"),
##'                    prop=c(1, 1, 0.5, 0.5),
##'                    yield=c(50, 40, 25, 22)))
##' RY(dat)
##' @export
RY <- function(dat, colIDstand = "ID", colIDfocal = "focal", colProp = "prop",
               colY = "yield", avgIfReps = FALSE) {
  stopifnot(
    is.data.frame(dat),
    all(c(colIDstand, colIDfocal, colProp, colY) %in% colnames(dat)),
    is.numeric(dat[[colProp]]),
    is.numeric(dat[[colY]]),
    all(!is.na(dat[[colProp]])),
    all(!is.na(dat[[colY]]))
  )

  out <- data.frame(dat,
    RY = NA
  )

  ## monovarietal cultures
  isMono <- (dat[[colProp]] == 1)
  stopifnot(
    sum(isMono) > 0,
    sum(isMono) < nrow(dat)
  )
  out$RY[which(isMono)] <- 1

  ## mixtures
  for (i in which(!isMono)) {
    mixYield <- dat[[colY]][i]
    focalID <- as.character(dat[[colIDfocal]][i])
    idx <- which(dat[[colIDstand]] == focalID &
      !is.na(dat[[colY]]))
    if (length(idx) > 1) {
      msg <- paste0("more than one monovarietal culture value for ", focalID)
      if (avgIfReps) {
        warning(msg, immediate. = TRUE)
      } else {
        stop(msg)
      }
    } else if (length(idx) == 0) {
      msg <- paste0("no monovarietal culture value for ", focalID)
      warning(msg, immediate. = TRUE)
      next
    }
    monoYield <- mean(dat[[colY]][idx], na.rm = TRUE)
    out$RY[i] <- mixYield / monoYield
  }

  return(out)
}

##' Relative yield total (RYT)
##'
##' Computes the "relative yield total" (RYT) of a mixture from \href{http://edepot.wur.nl/197683}{de Wit and Den Bergh (1965)} as the sum of all relative yields (RY) of the mixture.
##' It corresponds to equation 26 (based on 25a and not 25b!) in Weigelt and Jolliffe (2003, \doi{10.1046/j.1365-2745.2003.00805.x}), and is equivalent to the index called "land equivalent ratio" (LER) corresponding to equation 28 of the same paper.
##' Missing values are forbidden.
##' @param dat data.frame with one yield value per row and with at least four columns: the identifier of each stand, the identifier of the focal genotype or species in each stand, the sowing (plant) proportion of the focal in each stand, and the focal yield in each stand (see the example below); a stand with a single proportion equals to 1 will be interpreted as a monovarietal culture ("pure stands") and a mixture otherwise; note that contrary to \code{\link{RYP}}, the proportions will not be used for anything else
##' @param colIDstand column name for the stand identifiers
##' @param colIDfocal column name for the focal identifiers
##' @param colProp column name for the proportions
##' @param colY column name for the yield values
##' @return input data.frame with an extra column named "RYT"
##' @seealso \code{\link{RYP}}, \code{\link{RY}}, \code{\link{RYM}}, \code{\link{LER}}
##' @examples
##' (dat <- data.frame(ID=c("geno1", "geno2", "mixg1g2", "mixg1g2"),
##'                    focal=c("geno1", "geno2", "geno1", "geno2"),
##'                    prop=c(1, 1, 0.5, 0.5),
##'                    yield=c(50, 40, 25, 22)))
##' RYT(dat)
##' @export
RYT <- function(dat, colIDstand = "ID", colIDfocal = "focal", colProp = "prop",
                colY = "yield") {
  out <- RY(
    dat = dat, colIDstand = colIDstand, colIDfocal = colIDfocal,
    colProp = colProp, colY = colY
  )

  out <- data.frame(out,
    RYT = NA
  )

  isMix <- (dat[[colProp]] < 1)
  mixIDs <- unique(as.character(dat[[colIDstand]][which(isMix)]))
  for (mixID in mixIDs) {
    idx <- which(dat[[colIDstand]] == mixID)
    stopifnot(length(idx) > 1)
    out$RYT[idx] <- sum(out$RY[idx])
  }

  return(out)
}

##' Relative yield of mixture (RYM)
##'
##' Computes the "relative yield of mixture" (RYM) as defined initially by Wilson (1988, \doi{10.2307/2403626}) for binary mixtures sowed in equal proportions, and extended by Williams and McCarthy (2001, \doi{10.1046/j.1440-1703.2001.00368.x}) to binary mixtures of unequal proportions.
##' A RYM above 1 means that the mixture yielded better than the sum of the pure-stand yields, each weighted by their proportion in the mixture.
##' It corresponds to equation 35 in Weigelt and Jolliffe (2003, \doi{10.1046/j.1365-2745.2003.00805.x}), as well as equation 25b (and not to 25a!) of the same paper where the "control" treatment is the expected mixture yield calculated based on the weighted component monovarietal yields.
##' Note that the RYM was unfortunately called "relative yield" (RY) by Reiss and Drinkwater (2018, \doi{10.1002/eap.1629}).
##' Missing values are forbidden.
##' @param dat data.frame with one yield value per row and with at least four columns: the identifier of each stand, the identifiers of the components of each stand (separated by a specific symbol in the case of mixtures), the sowing (plant) proportion of the components in each stand (separated by the same symbol in the case of mixtures, and the yield of each stand (see the example below); a stand with a single component should have its proportion equal to 1 and will be interpreted as a monovarietal culture ("pure stand"), a mixture otherwise (and its proportions should sum to 1)
##' @param colIDstand column name for the stand identifiers (they should be unique); if NULL, the code will attempt to use \code{"colIDcomps"}
##' @param colIDcomps column name for the identifiers of the stand components (separated by \code{sep})
##' @param colProps column name for the stand proportions (separated by \code{sep} and summing to 1); if NULL, the code will assume all components of a given mixture are equiprobable
##' @param colY column name for the yield values
##' @param sep separator
##' @param colOut column name for the output RYM
##' @return input data.frame with an extra column named from \code{"colOut"}
##' @seealso \code{\link{estimRYMRep}}, \code{\link{RYP}}, \code{\link{RY}}, \code{\link{RYT}}
##' @examples
##' (dat <- data.frame(ID=c("geno1", "geno2", "mixg1g2"),
##'                    comps=c("geno1", "geno2", "geno1-geno2"),
##'                    props=c("1", "1", "0.6-0.4"),
##'                    yield=c(50, 40, 47)))
##' RYM(dat)
##' @export
RYM <- function(dat, colIDstand = NULL, colIDcomps = "comps", colProps = NULL,
                colY = "yield", sep = "-", colOut = "RYM") {
  stopifnot(
    is.data.frame(dat),
    all(c(colIDcomps, colY) %in% colnames(dat)),
    is.numeric(dat[[colY]]),
    all(!is.na(dat[[colY]]))
  )
  if (!is.null(colIDstand)) {
    stopifnot(
      colIDstand %in% colnames(dat),
      all(table(dat[[colIDstand]]) == 1)
    )
  }
  if (!is.null(colProps)) {
    stopifnot(all(!is.na(dat[[colProps]])))
  }

  out <- data.frame(dat)
  out[[colOut]] <- NA

  comps <- strsplit(as.character(dat[[colIDcomps]]), sep, fixed = TRUE)
  if (!is.null(colIDstand)) {
    IDs <- as.character(dat[[colIDstand]])
  } else {
    IDs <- as.character(dat[[colIDcomps]])
    if (!all(table(IDs) == 1)) {
      msg <- paste0("'colIDstand' is missing and IDs made from 'colIDcomps' are not unique")
      stop(msg)
    }
  }
  names(comps) <- IDs

  isMix <- (sapply(comps, length) > 1)
  if (sum(isMix) == nrow(dat)) {
    msg <- "each row of 'dat' is a mixture, monovarietal cultures are missing"
    stop(msg)
  }
  idxMixes <- which(isMix)
  mixIDs <- IDs[which(isMix)]
  stopifnot(all(!duplicated(mixIDs)))

  if (!is.null(colProps)) {
    props <- lapply(
      strsplit(as.character(dat[[colProps]]), sep, fixed = TRUE),
      as.numeric
    )
    names(props) <- IDs
  } else {
    props <- lapply(comps, function(x) {
      rep(1 / length(x), length(x))
    })
    names(props) <- IDs
  }
  stopifnot(
    all(sapply(props[idxMixes], length) > 1),
    all(unlist(Map(
      function(x) {
        all.equal(x, 1, tolerance = 10^(-6))
      },
      sapply(props[idxMixes], sum)
    )))
  )

  for (mixID in mixIDs) {
    idxMix <- which(IDs == mixID)
    stopifnot(length(idxMix) == 1)
    mixYield <- dat[[colY]][idxMix]
    mixComps <- comps[[mixID]]
    idxMonos <- which(IDs %in% mixComps)
    stopifnot(length(idxMonos) > 1)
    monoYields <- setNames(dat[[colY]][idxMonos], IDs[idxMonos])
    stopifnot(
      length(monoYields) == length(mixComps),
      all(!is.na(monoYields))
    )
    propComps <- setNames(props[[mixID]], mixComps)
    stopifnot(
      length(propComps) == length(monoYields),
      all(sort(names(propComps)) == sort(names(monoYields)))
    )
    propComps <- propComps[names(monoYields)] # force the same order
    out[[colOut]][idxMix] <- mixYield / (propComps %*% monoYields)[1, 1]
  }

  return(out)
}

##' Land equivalent ratio (LER)
##'
##' Computes the land equivalent ratio (LER) of Willey and Osiru (1972, \doi{10.1017/S0021859600025909}) interpreted as the relative land area required under sole cropping to produce the same yield as under intercropping.
##' It corresponds to equation 28 in Weigelt and Jolliffe (2003, \doi{10.1046/j.1365-2745.2003.00805.x}), and is equivalent to the index called "relative yield total" (RYT) corresponding to equation 26 of the same paper.
##' Missing values are forbidden.
##' @param x data frame with species in rows (with species names as row names) and at least two columns, the first being the performance under sole-cropping and the second under inter-cropping
##' @return list with partial and total LERs
##' @seealso \code{\link{RYT}}, \code{\link{RY}}
##' @author Timothee Flutre
##' @examples
##' (dat <- data.frame(solecrop=c(5, 15),
##'                    intercrop=c(4, 9),
##'                    row.names=c("grain", "fruit")))
##' LER(dat)
##' @export
LER <- function(x) {
  stopifnot(
    !is.null(rownames(x)),
    ncol(x) >= 2,
    all(!is.na(x))
  )

  species <- rownames(x)
  out <- list(pLER = stats::setNames(x[, 2] / x[, 1], species))
  out$LER <- sum(out$pLER)

  return(out)
}

##' Change in contribution
##'
##' Computes the change in contribution (CC) proposed by Williams and McCarthy (2001, \doi{10.1046/j.1440-1703.2001.00368.x}) (this function implements the equation for mixtures with any number of components as written in appendix V of the paper).
##' @param dat data frame with one yield value per row and with at least four columns: the identifier of each stand, the identifier of the focal genotype or species in each stand, the sowing (plant) proportion of the focal in each stand, and the focal yield in each stand (see the example below); a stand with a single proportion equals to 1 will be interpreted as a monovarietal culture ("pure stand"), a mixture otherwise (and its proportions should sum to 1)
##' @param colIDstand column name for the stand identifiers
##' @param colIDfocal column name for the focal identifiers
##' @param colProp column name for the proportions
##' @param colY column name for the yield values
##' @return input data.frame with an extra column named "CC"
##' @author Timothee Flutre
##' @examples
##' (dat <- data.frame(ID=c("geno1", "geno2", "mixg1g2", "mixg1g2"),
##'                    focal=c("geno1", "geno2", "geno1", "geno2"),
##'                    prop=c(1, 1, 0.5, 0.5),
##'                    yield=c(50, 40, 25, 22)))
##' CC(dat)
##' @export
CC <- function(dat, colIDstand = "ID", colIDfocal = "focal", colProp = "prop",
               colY = "yield") {
  stopifnot(
    is.data.frame(dat),
    all(c(colIDstand, colIDfocal, colProp, colY) %in% colnames(dat)),
    is.numeric(dat[[colProp]]),
    is.numeric(dat[[colY]]),
    all(!is.na(dat[[colProp]])),
    all(!is.na(dat[[colY]]))
  )

  out <- data.frame(dat,
    CC = NA
  )

  ## monovarietal cultures
  isMono <- (dat[[colProp]] == 1)
  stopifnot(
    sum(isMono) > 0,
    sum(isMono) < nrow(dat)
  )

  ## mixtures
  isMix <- (!isMono)
  idxMixes <- which(isMix)
  mixSumProps <- tapply(
    dat[[colProp]][idxMixes], dat[[colIDstand]][idxMixes],
    sum
  )
  stopifnot(all(unlist(Map(function(x) {
    all.equal(x, 1)
  }, mixSumProps))))
  mixIDs <- unique(as.character(dat[[colIDstand]][idxMixes]))
  for (mixID in mixIDs) {
    idxMix <- which(dat[[colIDstand]] == mixID)
    focalIDs <- as.character(dat[[colIDfocal]][idxMix])
    names(idxMix) <- focalIDs
    mixYields <- setNames(dat[[colY]][idxMix], focalIDs)
    props <- setNames(dat[[colProp]][idxMix], focalIDs)
    idxMono <- which(dat[[colIDstand]] %in% focalIDs)
    names(idxMono) <- dat[[colIDstand]][idxMono]
    idxMono <- idxMono[focalIDs]
    monoYields <- setNames(dat[[colY]][idxMono], names(idxMono))
    monoYields <- monoYields[focalIDs]
    for (i in seq_along(focalIDs)) {
      focalID <- focalIDs[i]
      out$CC[idxMix[focalID]] <- (mixYields[focalID] / sum(mixYields[-i])) /
        ((props[focalID] * monoYields[focalID]) / (sum(props[-i] * monoYields[-i]))) - 1
    }
  }

  return(out)
}

##' Index from Reiss and Drinkwater (2018)
##'
##' Computes the index of Reiss and Drinkwater (2018, \doi{10.1002/eap.1629}) that the authors unfortunately called "relative yield" (RY) whereas it corresponds to the "relative yield of mixture" (RYM).
##' Missing values are forbidden.
##' @param mixYields named numeric vector (or 1D-array) of mixture performances
##' @param monoYields named numeric vector (or 1D-array) of performances of the corresponding pure stands
##' @param mix2pur named list which components are the proportions of all the components for each mixture
##' @return numeric matrix with three columns, performance in the mixed stands, mean of their pure stands and overyielding, and as many rows as mixed stands
##' @seealso \code{\link{RYP}}, \code{\link{RY}}, \code{\link{RYT}}, \code{\link{RYM}}
##' @examples
##' mixYields <- setNames(c(40, 50),
##'                     c("mix2", "mix1"))
##' monoYields <- setNames(c(70, 35, 20),
##'                     c("varA1", "varC2", "varB3"))
##' mix2pur <- list("mix1"=setNames(c(0.5, 0.5), c("varB3", "varC2")),
##'                 "mix2"=setNames(c(0.5, 0.2, 0.3), c("varC2", "varA1", "varB3")))
##' RY_RD18(mixYields, monoYields, mix2pur)
##' @author Timothee Flutre
##' @noRd
RY_RD18 <- function(mixYields, monoYields, mix2pur) {
  stopifnot(
    is.numeric(mixYields),
    !is.null(names(mixYields)),
    is.numeric(monoYields),
    !is.null(names(monoYields)),
    is.list(mix2pur),
    !is.null(names(mix2pur)),
    all(names(mixYields) %in% names(mix2pur)),
    all(sapply(mix2pur, sum) == 1)
  )

  out <- lapply(names(mixYields), function(mix) {
    rowVecPerPur <- monoYields[names(mix2pur[[mix]])]
    colVecProp <- matrix(mix2pur[[mix]], dimnames = list(names(mix2pur[[mix]])))
    stopifnot(all(names(rowVecPerPur) == names(colVecProp)))
    weightedPurs <- (rowVecPerPur %*% colVecProp)[1, 1]
    RYM <- mixYields[mix] / weightedPurs
    names(RYM) <- NULL
    c(mixYields[mix], weightedPurs, RYM)
  })
  out <- do.call(rbind, out)
  rownames(out) <- names(mixYields)
  colnames(out) <- c("mix", "weightedPurs", "RY_RD18")

  return(out)
}

##' Overyielding
##'
##' Computes an index sometimes called overyielding (OY) but also known as "relative yield of mixture".
##' It hence is advised to use the function \code{\link{RYM}} because OY < 1 in fact corresponds to "underyielding".
##' Missing values are ignored.
##' @param perfMix named numeric vector (or 1D-array) of mixture performances
##' @param perfPur named numeric vector (or 1D-array) of performances of the corresponding pure stands
##' @param mix2pur named list which components are the mixture composition
##' @return numeric matrix with three columns, performance in the mixed stands, mean of their pure stands and overyielding, and as many rows as mixed stands
##' @author Timothee Flutre
##' @examples
##' perfMix <- setNames(c(40, 50),
##'                     c("mix2", "mix1"))
##' perfPur <- setNames(c(70, 35, 20),
##'                     c("varA1", "varC2", "varB3"))
##' mix2pur <- list("mix1"=c("varB3", "varC2"),
##'                 "mix2"=c("varC2", "varA1", "varB3"))
##' plantmix:::overyielding(perfMix, perfPur, mix2pur)
##' @noRd
overyielding <- function(perfMix, perfPur, mix2pur) {
  stopifnot(
    is.numeric(perfMix),
    !is.null(names(perfMix)),
    is.numeric(perfPur),
    !is.null(names(perfPur)),
    is.list(mix2pur),
    !is.null(names(mix2pur)),
    all(names(perfMix) %in% names(mix2pur))
  )

  out <- lapply(names(perfMix), function(mix) {
    if (any(is.na(perfPur[mix2pur[[mix]]]))) {
      msg <- paste0(
        "missing performance in pure stand for genotype(s)",
        " in mixture '", mix, "'"
      )
      warning(msg)
    }
    meanPur <- mean(perfPur[mix2pur[[mix]]], na.rm = TRUE)
    OY <- perfMix[mix] / meanPur
    names(OY) <- NULL
    c(perfMix[mix], meanPur, OY)
  })
  out <- do.call(rbind, out)
  rownames(out) <- names(perfMix)
  colnames(out) <- c("mix", "meanPur", "OY")

  return(out)
}


##' Relative interaction index (RII)
##'
##' Computes the relative interaction index (RII) as defined in Armas et al (2004, \doi{10.1890/03-0650}): RII_i = (Y_i(j) - Y_i) / (Y_i(j) + Y_i) where Y_i(j) is the yield per plant of i when mixed with j, and Y_i is the yield per plant of i when growned in isolation.
##' Missing values are forbidden.
##' @param dat data.frame with one yield value per row and with at least four columns: the identifier of each stand, the identifier of the focal genotype or species in each stand, the sowing (plant) proportion of the focal in each stand, and the focal yield in each stand (see the example below); a stand with a single proportion equals to 1 will be interpreted as a monovarietal culture ("pure stand"), a mixture otherwise (and its proportions should sum to 1)
##' @param colIDstand column name for the stand identifiers
##' @param colIDfocal column name for the focal identifiers
##' @param colProp column name for the proportions
##' @param colY column name for the yield values
##' @param avgIfReps if TRUE, replicates of monovarietals will be averaged; if FALSE, the presence of replicates will return an error
##' @return input data.frame with an extra column named "RII"
##' @author Timothee Flutre
##' @seealso \code{\link{RIInet}}, \code{\link{RY}}, \code{\link{RYT}}, \code{\link{RYP}}, \code{\link{RYM}}
##' @examples
##' sow_density <- 200
##' (dat <- data.frame(ID=c("geno1", "geno2", "mixg1g2", "mixg1g2"),
##'                    focal=c("geno1", "geno2", "geno1", "geno2"),
##'                    prop=c(1, 1, 0.5, 0.5),
##'                    yield_qt_ha=c(50, 40, 25, 22)))
##' dat$yield_g_m2 <- (dat$yield_qt_ha / 10^4) * 10^6
##' dat$yield_g_plant <- dat$yield_g_m2 / (sow_density * dat$prop)
##' RII(dat, colY = "yield_g_plant")
##' @export
RII <- function(dat, colIDstand = "ID", colIDfocal = "focal", colProp = "prop",
                colY = "yield", avgIfReps = FALSE) {
  stopifnot(
    is.data.frame(dat),
    all(c(colIDstand, colIDfocal, colProp, colY) %in% colnames(dat)),
    is.numeric(dat[[colProp]]),
    is.numeric(dat[[colY]]),
    all(!is.na(dat[[colProp]])),
    all(!is.na(dat[[colY]]))
  )

  out <- data.frame(dat,
    RII = NA
  )

  ## monovarietal cultures
  isMono <- (dat[[colProp]] == 1)
  stopifnot(
    sum(isMono) > 0,
    sum(isMono) < nrow(dat)
  )
  out$RII[which(isMono)] <- 1

  ## mixtures
  for (i in which(!isMono)) {
    mixYield <- dat[[colY]][i]
    focalID <- as.character(dat[[colIDfocal]][i])
    idx <- which(dat[[colIDstand]] == focalID &
      !is.na(dat[[colY]]))
    if (length(idx) > 1) {
      msg <- paste0("more than one monovarietal culture value for ", focalID)
      if (avgIfReps) {
        warning(msg, immediate. = TRUE)
      } else {
        stop(msg)
      }
    } else if (length(idx) == 0) {
      msg <- paste0("no monovarietal culture value for ", focalID)
      warning(msg, immediate. = TRUE)
      next
    }
    monoYield <- mean(dat[[colY]][idx], na.rm = TRUE)
    out$RII[i] <- (mixYield - monoYield) / (mixYield + monoYield)
  }

  return(out)
}

##' Net Relative interaction index (RIInet)
##'
##' Computes the net relative interaction index (RII) as defined in \href{https://elifesciences.org/articles/77577}{Stefan et al (2022)}: RIInet = sum_i prop_i RII_i.
##' Missing values are forbidden.
##' @param dat data.frame with one yield value per row and with at least four columns: the identifier of each stand, the identifier of the focal genotype or species in each stand, the sowing (plant) proportion of the focal in each stand, and the relative interaction index of the focal.
##' @param colIDstand column name for the stand identifiers
##' @param colIDfocal column name for the focal identifiers
##' @param colProp column name for the proportions
##' @param colRII column name for the RII values
##' @return data.frame with columns \code{colIDstand} and "RIInet"
##' @author Timothee Flutre
##' @seealso \code{\link{RII}}
##' @examples
##' sow_density <- 200
##' (dat <- data.frame(ID=c("geno1", "geno2", "mixg1g2", "mixg1g2"),
##'                    focal=c("geno1", "geno2", "geno1", "geno2"),
##'                    prop=c(1, 1, 0.5, 0.5),
##'                    yield_qt_ha=c(50, 40, 25, 22)))
##' dat$yield_g_m2 <- (dat$yield_qt_ha / 10^4) * 10^6
##' dat$yield_g_plant <- dat$yield_g_m2 / (sow_density * dat$prop)
##' dat <- RII(dat, colY = "yield_g_plant")
##' RIInet(dat)
##' @export
RIInet <- function(dat, colIDstand = "ID", colIDfocal = "focal", colProp = "prop",
                   colRII = "RII") {
  stopifnot(
    is.data.frame(dat),
    all(c(colIDstand, colIDfocal, colProp, colRII) %in% colnames(dat)),
    is.numeric(dat[[colProp]]),
    is.numeric(dat[[colRII]]),
    all(!is.na(dat[[colProp]])),
    all(!is.na(dat[[colRII]]))
  )

  ## monovarietal cultures
  isMono <- (dat[[colProp]] == 1)
  if (any(isMono)) {
    dat <- droplevels(dat[!isMono, ])
  }

  ## mixtures
  mixIDs <- NA
  if (is.factor(dat[[colIDstand]])) {
    mixIDs <- levels(dat[[colIDstand]])
  } else {
    mixIDs <- unique(sort(as.character(dat[[colIDstand]])))
  }

  out <- list()
  out[[colIDstand]] <- mixIDs
  out[["RIInet"]] <- rep(NA, length(mixIDs))
  out <- as.data.frame(out)

  out$RIInet <- Map(function(mixID) {
    idx <- grep(paste0("^", mixID, "$"), dat[[colIDstand]])
    sum(dat[idx, colProp] * dat[idx, colRII])
  }, mixIDs)
  out$RIInet <- unlist(out$RIInet)

  return(out)
}

## ---------------------------------------------------------------------------
## functions written during the internship of N. Vazeux-Blumental
## TODO: check if they correspond to indices in Weigelt and Jolliffe (2003)

##' Component ratio
##'
##' Computes the ratio of the output (e.g., yield) of each component in the mixture versus in pure stands.
##' @param yieldPure named numeric vector of performances of the corresponding pure stands
##' @param yieldMix named numeric vector of mixture performances
##' @param order_mix order of the mixture, (e.g., 2 for a binary mixture)
##' @return numeric vector
##' @examples
##' yieldPure <- c("var37"=12, "var12"=16)
##' yieldMix <- c("var12"=15, "var37"=14)
##' ratioComp(yieldPure, yieldMix)
##' @seealso \code{\link{ratioProd}}, \code{\link{ratioAggr}}, \code{\link{ratioRelAggr}}
##' @author Noa Vazeux-Blumental, Timothee Flutre
##' @noRd
ratioComp <- function(yieldPure, yieldMix, order_mix = NULL) {
  stopifnot(
    is.vector(yieldPure),
    is.numeric(yieldPure),
    is.vector(yieldMix),
    is.numeric(yieldMix)
  )
  if (all(!is.null(names(yieldPure)), !is.null(yieldMix))) {
    stopifnot(all(sort(names(yieldMix)) == sort(names(yieldPure))))
    yieldMix <- yieldMix[names(yieldPure)]
  }

  out <- c()

  order_mix <- length(yieldMix)
  for (i in 1:length(yieldPure)) {
    out[i] <- c(yieldMix[i] * order_mix / yieldPure[i])
  }
  if (all(!is.null(names(yieldPure)), !is.null(names(yieldMix)))) {
    names(out) <- names(yieldPure)
  }

  return(out)
}

##' Productivity ratio
##'
##' Computes the productivity ratio of the mixture, similar to the overyielding.
##' @param ratioComp_vect output from \code{\link{ratioComp}}
##' @return numeric
##' @seealso \code{\link{ratioComp}}, \code{\link{ratioAggr}}, \code{\link{ratioRelAggr}}
##' @author Noa Vazeux-Blumental, Timothee Flutre
##' @noRd
ratioProd <- function(ratioComp_vect) {
  stopifnot(
    is.vector(ratioComp_vect),
    is.numeric(ratioComp_vect)
  )
  order_mix <- length(ratioComp_vect)
  out <- sum(ratioComp_vect) / order_mix
  return(out)
}

##' Aggressivity ratio
##'
##' Computes the aggressivity ratio (only for binary mixtures).
##' @param ratioComp_vect output from \code{\link{ratioComp}}
##' @return numeric vector
##' @seealso \code{\link{ratioComp}}, \code{\link{ratioProd}}, \code{\link{ratioRelAggr}}
##' @author Noa Vazeux-Blumental, Timothee Flutre
##' @noRd
ratioAggr <- function(ratioComp_vect) {
  stopifnot(
    is.vector(ratioComp_vect),
    is.numeric(ratioComp_vect),
    length(ratioComp_vect) == 2
  )
  order_mix <- 2
  out <- rep(NA, order_mix)
  out[1] <- ratioComp_vect[1] - ratioComp_vect[2]
  out[2] <- ratioComp_vect[2] - ratioComp_vect[1]
  return(out)
}

##' Relative aggressivity ratio
##'
##' Computes the relative aggressivity ratio (only for binary mixtures).
##' @param ratioComp_vect output from \code{\link{ratioComp}}
##' @return numeric vector
##' @seealso \code{\link{ratioComp}}, \code{\link{ratioProd}}, \code{\link{ratioAggr}}
##' @author Noa Vazeux-Blumental, Timothee Flutre
##' @noRd
ratioRelAggr <- function(ratioComp_vect) {
  stopifnot(
    is.vector(ratioComp_vect),
    is.numeric(ratioComp_vect),
    length(ratioComp_vect) == 2
  )
  order_mix <- 2
  out <- rep(NA, order_mix)
  out[1] <- ratioComp_vect[1] / ratioComp_vect[2]
  out[2] <- ratioComp_vect[2] / ratioComp_vect[1]
  return(out)
}
