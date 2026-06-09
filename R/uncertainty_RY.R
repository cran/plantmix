## Functions implementing methods to quantify the uncertainty around RY.

##' RY with "rep" effect
##'
##' Computes the relative yields (RYs) per replicate (e.g., blocks), then fits a linear model correcting for this rep effect, and quantifies the uncertainty around it.
##' @param data data frame
##' @param colY column name specifying the response variable
##' @param colIDstand column name specifying the stand identifiers
##' @param colIDfocal column name for the focal identifiers
##' @param colRep column name specifying the replicate
##' @param mix2genos named list with one component per mixture, each component being a vector of genotype names; see \code{\link{getMixtureList}}
##' @param plotDiag if TRUE, diagnostic plots
##' @param conf_level confidence level for the uncertainty intervals
##' @return data frame with the inference summaries of RY per mixture
##' @author Timothee Flutre [aut], Arnaud Gauffreteau [ctb]
##' @seealso \code{\link{estimRYMRep}}
##' @examples
##' (dat <- data.frame(ID=c("geno1", "geno1",
##'                         "geno2", "geno2",
##'                         "mixg1g2", "mixg1g2", "mixg1g2", "mixg1g2"),
##'                    focal=c("geno1", "geno1",
##'                            "geno2", "geno2",
##'                            "geno1", "geno2",
##'                            "geno1", "geno2"),
##'                    prop=c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
##'                    block=c("A","B","A","B","A","A","B","B"),
##'                    yield=c(51,45, 39,43, 25,22,26,25)))
##' estimRYRep(dat, "yield", "ID", "focal", "block", list("mixg1g2"=c("geno1","geno2")))
##' @export
estimRYRep <- function(data, colY, colIDstand, colIDfocal, colRep, mix2genos,
                       conf_level = 0.95, plotDiag = FALSE) {
  if (FALSE) { # debug
    data <- dat
    colY <- "yield"
    colIDstand <- "ID"
    colIDfocal <- "focal"
    colRep <- "block"
    mix2genos <- list("mixg1g2" = c("geno1", "geno2"))
    conf_level <- 0.95
    plotDiag <- FALSE
  }
  stopifnot(
    requireNamespace("emmeans", quietly = TRUE),
    is.data.frame(data),
    all(c(colY, colIDstand, colIDfocal, colRep) %in% names(data)),
    is.list(mix2genos),
    !is.null(names(mix2genos)),
    conf_level >= 0 & conf_level <= 1,
    is.logical(plotDiag)
  )

  ## compute the RY of each mix per replicate (block)
  if (!is.factor(data[[colRep]])) {
    data[[colRep]] <- factor(data[[colRep]])
  }
  repData <- lapply(levels(data[[colRep]]), function(repLev) {
    subData <- droplevels(data[data[[colRep]] == repLev, ])
    n <- sum(sapply(mix2genos, length))
    repRYs <- list()
    repRYs[[colIDstand]] <- rep(names(mix2genos), sapply(mix2genos, length))
    repRYs[[colIDfocal]] <- do.call(c, mix2genos)
    repRYs[[colRep]] <- rep(repLev, n)
    repRYs$RY <- rep(NA, n)
    repRYs <- as.data.frame(do.call(cbind, repRYs))
    repRYs$RY <- as.numeric(repRYs$RY)
    for (i in 1:nrow(repRYs)) {
      mix <- repRYs[[colIDstand]][i]
      comps <- mix2genos[[mix]]
      for (focal in comps) {
        idx <- which(subData[[colIDstand]] == mix &
          subData[[colIDfocal]] == focal)
        if (length(idx) == 0) {
          next
        }
        if (length(idx) > 1) {
          msg <- paste0(
            length(idx), " data will be averaged for ", mix,
            " in replicate ", repLev
          )
          warning(msg, immediate. = TRUE)
        }
        mixY <- mean(subData[idx, colY])
        idx <- which(subData[[colIDstand]] == focal &
          subData[[colIDfocal]] == focal)
        if (length(idx) == 0) {
          next
        }
        if (length(idx) > 1) {
          msg <- paste0(
            length(idx), " data will be averaged for ", focal,
            " in replicate ", repLev
          )
          warning(msg, immediate. = TRUE)
        }
        monoY <- mean(subData[idx, colY])
        repRYs$RY[i] <- mixY / monoY
      }
    }
    repRYs
  })
  repData <- do.call(rbind, repData)
  repData[[colIDstand]] <- factor(repData[[colIDstand]])
  repData[[colIDfocal]] <- factor(repData[[colIDfocal]])
  repData[[colRep]] <- factor(repData[[colRep]])

  ## Fit a linear model correcting for the replicate effect
  repData$focal_stand <- paste0(repData[[colIDfocal]], " in ", repData[[colIDstand]])
  form <- as.formula(paste0("RY ~ ", colRep, " + focal_stand"))
  contr <- list()
  contr[[colRep]] <- "contr.sum"
  contr[["focal_stand"]] <- "contr.sum"
  fit <- lm(form, data = repData, contrasts = contr)
  if (plotDiag) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(2, 2))
    plot(fit, which = 1:4)
  }
  fitMM <- emmeans::emmeans(fit, specs = "focal_stand", level = conf_level)
  ## the next cmd is commented because the null may change depending on the mix
  ## fitTest <- as.data.frame(emmeans::test(fitMM, null = 1))
  fitMM <- as.data.frame(fitMM)
  fitMM[[colIDfocal]] <- sapply(
    strsplit(as.character(fitMM$focal_stand), " in "),
    `[`, 1
  )
  fitMM[[colIDstand]] <- sapply(
    strsplit(as.character(fitMM$focal_stand), " in "),
    `[`, 2
  )
  fitMM$estim <- fitMM$emmean
  fitMM$se <- fitMM$SE
  fitMM$cil <- fitMM$lower.CL
  fitMM$ciu <- fitMM$upper.CL
  fitMM$pv <- NA
  fitMM <- fitMM[, c(
    colIDstand, colIDfocal, "estim", "se", # "df",
    "cil", "ciu"
  )] # , "pv")]

  return(fitMM)
}
