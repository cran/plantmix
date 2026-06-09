## Functions implementing methods to quantify the uncertainty around RYM.

##' RYM with "rep" effect
##'
##' Computes the relative yield of mixtures (RYMs) per replicate (e.g., blocks), then fits a linear model correcting for this rep effect, and quantifies the uncertainty around it.
##' @param data data frame
##' @param colY column name specifying the response variable
##' @param colID column name specifying the stand identifiers
##' @param colRep column name specifying the replicate
##' @param mix2genos named list with one component pe rmixture, each component being a vector of genotype names
##' @param conf_level confidence level for the uncertainty intervals
##' @param plotDiag if TRUE, diagnostics will be plotted
##' @param title title for the diagnostics plots
##' @return data frame with the inference summaries of RYM per mixture
##' @author Timothee Flutre [aut], Arnaud Gauffreteau [ctb]
##' @seealso \code{\link{estimRYRep}}
##' @examples
##' ## generate fake data
##' set.seed(1234)
##' nbGenos <- 25
##' genos <- sprintf("g%02i", 1:nbGenos)
##' pairs <- t(combn(x=genos, m=2))
##' mixIDs <- sample(paste(pairs[,1], pairs[,2], sep="_"), size=75)
##' nbBlocks <- 3
##' blocks <- LETTERS[1:nbBlocks]
##' dat <- do.call(rbind, lapply(blocks, function(block){
##'   data.frame(ID=c(genos, mixIDs),
##'              block=block,
##'              stringsAsFactors=TRUE)
##' }))
##' listContr <- list(block="contr.sum")
##' X <- model.matrix(~ 1 + block, data=dat, contrasts=listContr)
##' Z_GMA <- mkZGMA(dat, "ID", "_")
##' Z_SMA <- mkZSMA(dat, "ID", "_", inc_SMA_ii="only_pur")
##' truth <- list("intercept"=100, "var_GMA"=10, "var_SMA"=4, "var_error"=1)
##' truth[["blockEffs"]] <- rnorm(n=nbBlocks - 1, mean=0, sd=2)
##' truth[["GMAs"]] <- rnorm(n=nbGenos, mean=0, sd=sqrt(truth$var_GMA))
##' truth[["SMAs"]] <- rnorm(n=nlevels(dat$ID), mean=0, sd=sqrt(truth$var_SMA))
##' truth[["errors"]] <- rnorm(n=nrow(dat), mean=0, sd=sqrt(truth$var_error))
##' y <- X %*% c(truth$intercept, truth$blockEffs) +
##'   Z_GMA %*% truth$GMAs +
##'   Z_SMA %*% truth$SMAs +
##'   truth$errors
##' dat$yield <- y[,1]
##' hist(dat$yield, breaks=20, las=1, main="Simulated data")
##' boxplot(yield ~ block, data=dat, las=1, main="Simulated data")
##'
##' ## compute the average yield per ID over blocks, and then compute the RYMs:
##' avg <- tapply(dat$yield, dat$ID, mean)
##' avg <- data.frame(ID=names(avg), yield=avg, stringsAsFactors=TRUE)
##' avg <- RYM(avg, colIDcomps="ID", colY="yield", sep="_")
##'
##' ## compute the RYMs per block, and then correct for the "block" effect:
##' mix2genos <- strsplit(levels(dat$ID), "_")
##' names(mix2genos) <- levels(dat$ID)
##' out <- estimRYMRep(dat, "yield", "ID", "block", mix2genos)
##' head(out)
##'
##' ## compare both approaches:
##' op <- par(mfrow=c(1,2))
##' hist(avg$RYM, las=1, main="Estimated RYM by averaging over blocks", xlab="RYM")
##' abline(v=1, lwd=2); abline(v=mean(avg$RYM, na.rm=TRUE), col="red", lwd=2)
##' hist(out$estim, las=1, main="Estimated RYM after correcting the 'block' effect", xlab="RYM")
##' abline(v=1, lwd=2); abline(v=mean(out$estim), col="red", lwd=2)
##' par(op)
##'
##' out$avgRYM <- avg[rownames(out), "RYM"]
##' plot(out$avgRYM, out$estim, las=1, type="n",
##'      xlab="RYM averaged over blocks",
##'      ylab="RYM after correcting the 'block' effect",
##'      main="Estimated RYMs")
##' idx <- which(out$pv > 0.05)
##' points(out$avgRYM[idx], out$estim[idx], pch=19, col="black")
##' idx <- which(out$pv <= 0.05)
##' points(out$avgRYM[idx], out$estim[idx], pch=19, col="red")
##' abline(a=0, b=1, h=1, v=1, lty=2)
##' segments(x0=as.numeric(out$avgRYM), y0=out$cil,
##'          x1=as.numeric(out$avgRYM), y1=out$ciu,
##'          col="grey", lty=2, lwd=2)
##' @export
estimRYMRep <- function(data, colY, colID, colRep, mix2genos, conf_level = 0.95,
                        plotDiag = FALSE, title = "") {
  stopifnot(
    requireNamespace("emmeans", quietly = TRUE),
    is.data.frame(data),
    all(c(colY, colID, colRep) %in% names(data)),
    is.list(mix2genos),
    !is.null(names(mix2genos)),
    conf_level >= 0 & conf_level <= 1,
    is.logical(plotDiag)
  )

  out <- data.frame(
    mix = names(mix2genos),
    comps = sapply(mix2genos, paste0, collapse = "-"),
    size = sapply(mix2genos, length),
    estim = NA,
    se = NA,
    cil = NA,
    ciu = NA,
    pv = NA
  )

  ## compute the RYM of each mix per replicate (block)
  repData <- lapply(levels(data[[colRep]]), function(repLev) {
    repRYMs <- data.frame(RYM = rep(NA, length(mix2genos)))
    repRYMs[[colRep]] <- repLev
    repRYMs[[colID]] <- names(mix2genos)
    subData <- droplevels(data[data[[colRep]] == repLev, ])
    for (i in 1:nrow(repRYMs)) {
      mix <- repRYMs[[colID]][i]
      idxMix <- which(subData[[colID]] == mix)
      respMix <- subData[[colY]][idxMix]
      respPurs <- rep(NA, length(mix2genos[[mix]]))
      for (j in seq_along(respPurs)) {
        pur <- mix2genos[[mix]][j]
        idxPur <- which(subData[[colID]] == pur)
        respPurs[j] <- subData[[colY]][idxPur]
      }
      repRYMs$RYM[i] <- respMix / mean(respPurs) # do NOT ignore NA in mean()
    }
    repRYMs
  })
  repData <- do.call(rbind, repData)
  repData[[colRep]] <- factor(repData[[colRep]])
  repData[[colID]] <- factor(repData[[colID]])

  ## Fit a linear model correcting for the replicate effect
  form <- as.formula(paste0("RYM ~ ", colRep, " + ", colID))
  contr <- list()
  contr[[colRep]] <- "contr.sum"
  contr[[colID]] <- "contr.sum"
  fit <- lm(form, data = repData, contrasts = contr)
  if (plotDiag) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(2, 2))
    plot(fit, which = 1:4, main = title)
  }
  fitMM <- emmeans::emmeans(fit, specs = colID, level = conf_level)
  fitTest <- as.data.frame(emmeans::test(fitMM, null = 1)) # H0: "RYM = 1"
  fitMM <- as.data.frame(fitMM)
  idxFit <- match(out$mix, fitMM[[colID]])
  idxOut <- which(!is.na(idxFit))
  idxFit <- idxFit[!is.na(idxFit)]
  stopifnot(all(as.character(out$mix[idxOut]) == as.character(fitMM[[colID]][idxFit])))
  out$estim[idxOut] <- fitMM$emmean[idxFit]
  out$se[idxOut] <- fitMM$SE[idxFit]
  out$cil[idxOut] <- fitMM$lower.CL[idxFit]
  out$ciu[idxOut] <- fitMM$upper.CL[idxFit]

  idxFit <- match(out$mix, fitTest[[colID]])
  idxOut <- which(!is.na(idxFit))
  idxFit <- idxFit[!is.na(idxFit)]
  stopifnot(all(as.character(out$mix[idxOut]) == as.character(fitTest[[colID]][idxFit])))
  ## out$estim[idxOut] <- fitTest$emmean[idxFit]
  ## out$se[idxOut] <- fitTest$SE[idxFit]
  out$pv[idxOut] <- fitTest$p.value[idxFit]

  return(out)
}
