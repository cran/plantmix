## functions useful to plan the sowing of a field trial

##' Sowing
##'
##' Computes the number of seeds of a given genotype to be sowned as pure stand based on a given sowing area and density.
##' @param sowingArea area to be sowned in square meters
##' @param sowingDensity number of seeds per square meter
##' @return number of seeds (in thousands)
##' @author Timothee Flutre
##' @seealso \code{\link{mixSowingWeight}}
##' @examples
##' sowingArea <- 9.52
##' sowingDensity <- 160
##' nbSeedsToSownInPure(sowingArea, sowingDensity)
##' sowingArea <- 8.66
##' sowingDensity <- 200
##' nbSeedsToSownInPure(sowingArea, sowingDensity)
##' @export
nbSeedsToSownInPure <- function(sowingArea = 8.4, sowingDensity = 160) {
  stopifnot(
    sowingDensity > 0,
    sowingArea > 0
  )
  out <- NA

  ## | area (m2)  |  nb seeds       |  weight in pure (g)     |
  ## |------------|-----------------|-------------------------|
  ## | 1          |  sowingDensity  |  sowingWeightOneSqMeter |
  ## |            | 1000            |  tkw                    |
  ## | sowingArea |                 |  sowingWeight
  tkw <- 45 # arbitrary because will be multiplied and then divided
  sowingWeightOneSqMeter <- (tkw * sowingDensity) / 1000
  sowingWeight <- (sowingArea * sowingWeightOneSqMeter) / 1

  out <- sowingWeight / tkw

  return(out)
}

##' Sowing
##'
##' Computes the seed weight per genotype for each stand, assuming equal sowing proportions.
##' @param pureSowingWeights vector which names are genotypes and values are sowing weights in pure stands
##' @param stands vector which names are stand identifiers and values are genotype(s), a single one if it is a pure stand, several ones if it is a mixed stand (separated by a dash "-", e.g., "var1-var8")
##' @return list which components are stands and values are named vectors with sowing weight per genotype for each stand
##' @author Timothee Flutre
##' @seealso \code{\link{nbSeedsToSownInPure}}
##' @examples
##' sowingArea <- 9.52
##' sowingDensity <- 160
##' (tmp <- nbSeedsToSownInPure(sowingArea, sowingDensity))
##' TKWs <- c("var1"=38.605, "var2"=40.051, "var6"=36.251, "var8"=33.368)
##' pureSowingWeights <- TKWs * tmp
##' stands <- c("var1"="var1", "mix8"="var1-var2", "mix34"="var1-var8",
##'             "mix50"="var1-var6-var8")
##' mixSowingWeight(pureSowingWeights, stands)
##' @export
mixSowingWeight <- function(pureSowingWeights, stands) {
  stopifnot(
    is.vector(pureSowingWeights),
    !is.null(names(pureSowingWeights)),
    is.vector(stands),
    !is.null(names(stands))
  )

  out <- list()
  for (stand in names(stands)) {
    genos <- strsplit(stands[[stand]], "-")[[1]]
    if (any(!genos %in% names(pureSowingWeights))) {
      msg <- paste0(
        "missing sowing weight in pure stand for genotype(s)",
        " of stand '", stand, "'"
      )
      warning(msg)
      next
    }
    standOrder <- length(genos) # 1 if pure stand, higher otherwise
    sowingProp <- 1 / standOrder
    out[[stand]] <- stats::setNames(rep(NA, standOrder), genos)
    for (geno in genos) {
      out[[stand]][geno] <- sowingProp * pureSowingWeights[geno]
    }
  }

  return(out)
}
