## Functions used in both optim_designs_var_mix.R and optim_designs_crop_mix.R

##' @noRd
convDiallel2Combs <- function(diallel) {
  ind <- which(upper.tri(diallel, diag = TRUE), arr.ind = TRUE)
  nn <- dimnames(diallel)
  out <- data.frame(
    row = nn[[1]][ind[, 1]],
    col = nn[[2]][ind[, 2]],
    val = diallel[ind],
    stringsAsFactors = FALSE
  )
  out <- out[out$val == 1, ]
  out$val <- NULL
  colnames(out) <- c("comp1", "comp2")
  return(out)
}

##' Plot diallel
##'
##' Plots a diallel matrix.
##' @param diallel matrix
##' @param main title
##' @return nothing
##' @author Timothee Flutre
##' @examples
##' nbGenos <- 25
##' levGenos <- sprintf(fmt=paste0("geno%0", floor(log10(nbGenos))+1, "i"),
##'                     1:nbGenos)
##' nbMixes <- 75 # only binary and balanced
##' design <- getDesignBinaryVarMix(levGenos, nbMixes, seed=12345)
##' plotDiallel(design$diallel)
##' @export
plotDiallel <- function(diallel, main = NULL) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(1, 4, 6, 1) + 0.1)
  graphics::image(t(diallel)[, nrow(diallel):1],
    axes = FALSE,
    col = c("white", "black")
  )
  if (is.null(main)) {
    combs <- convDiallel2Combs(diallel)
    idxS <- getShannonIndex(getGenoPropsFromCombs(combs))
    main <- ""
    if (nrow(diallel) == ncol(diallel)) {
      main <- paste0(nrow(diallel), " genotypes")
    } else {
      main <- paste0(nrow(diallel), " x ", ncol(diallel), " genotypes")
    }
    main <- paste0(
      main, "; ", sum(diallel), " binary mixtures",
      "; entropy = ", round(idxS, 4)
    )
  }
  title(main = main, line = 4)
  axis(side = 2, at = seq(0, 1, length.out = nrow(diallel)), labels = rev(rownames(diallel)), las = 1)
  axis(side = 3, at = seq(0, 1, length.out = ncol(diallel)), labels = colnames(diallel))
}
