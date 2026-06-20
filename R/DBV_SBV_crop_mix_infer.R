##' @noRd
fitDBVSBVinter_prepIn_listData <- function(listY, listX, listZ, listVCov) {
  ## small reformatting
  if (!is(listY$Y_IC, "matrix")) {
    listY$Y_IC <- as.matrix(listY$Y_IC)
  }
  if ("y_SC_f" %in% names(listY)) {
    if (is(listY$y_SC_f, "matrix")) {
      if (ncol(listY$y_SC_f) == 1) {
        listY$y_SC_f <- as.vector(listY$y_SC_f)
      }
    }
  }
  if ("y_SC_t" %in% names(listY)) {
    if (is(listY$y_SC_t, "matrix")) {
      if (ncol(listY$y_SC_t) == 1) {
        listY$y_SC_t <- as.vector(listY$y_SC_t)
      }
    }
  }

  listData <- list()
  listData$model <- "DBV_SBV_crop_mix"

  ## TMB does not allow DATA or PARAMETER in if statements,
  ## hence all variables must always be passed on;
  ## when you don't want to use one of them in the inference,
  ## set it to 0.0 but pass it on

  ## response variable(s)
  listData$Y_IC <- listY$Y_IC
  if ("y_SC_f" %in% names(listY)) {
    listData$y_SC_f <- listY$y_SC_f
  } else {
    listData$y_SC_f <- c(0.0)
  }
  if ("y_SC_t" %in% names(listY)) {
    listData$y_SC_t <- listY$y_SC_t
  } else {
    listData$y_SC_t <- c(0.0)
  }

  ## design matrices of fixed effects
  listData$X_IC <- sparseMat4Tmb(listX$X_IC)
  if ("X_SC_f" %in% names(listX)) {
    listData$X_SC_f <- sparseMat4Tmb(listX$X_SC_f)
  } else {
    listData$X_SC_f <- sparseMat4Tmb(matrix(0.0, 2, 2))
  }
  if ("X_SC_t" %in% names(listX)) {
    listData$X_SC_t <- sparseMat4Tmb(listX$X_SC_t)
  } else {
    listData$X_SC_t <- sparseMat4Tmb(matrix(0.0, 2, 2))
  }

  ## design matrices of random effects
  listData$Z_DS_f <- sparseMat4Tmb(listZ$Z_DS_f)
  if ("Z_D_f" %in% names(listZ)) {
    listData$Z_D_f <- sparseMat4Tmb(listZ$Z_D_f)
  } else {
    listData$Z_D_f <- sparseMat4Tmb(matrix(0.0, 2, 2))
  }
  if ("Z_DxS_f" %in% names(listZ)) {
    listData$Z_DxS_f <- sparseMat4Tmb(listZ$Z_DxS_f)
  } else {
    listData$Z_DxS_f <- sparseMat4Tmb(matrix(0.0, 2, 2))
  }

  ## kinship matrices
  listData$K <- listVCov$K[
    colnames(listData$Z_DS_f),
    colnames(listData$Z_DS_f)
  ]
  if ("Kmixpair" %in% names(listVCov)) {
    stopifnot(!"invKmixpair" %in% names(listVCov))
    listData$Kmixpair <- listVCov$Kmixpair[
      colnames(listData$Z_DxS_f),
      colnames(listData$Z_DxS_f)
    ]
  } else {
    listData$Kmixpair <- matrix(0.0, 2, 2)
  }
  if ("invKmixpair" %in% names(listVCov)) {
    stopifnot(!"Kmixpair" %in% names(listVCov))
    listData$invKmixpair <- listVCov$invKmixpair[
      colnames(listData$Z_DxS_f),
      colnames(listData$Z_DxS_f)
    ]
    if (!is(listData$invKmixpair, "sparseMatrix")) {
      listData$invKmixpair <- sparseMat4Tmb(listData$invKmixpair)
    }
  } else {
    listData$invKmixpair <- sparseMat4Tmb(matrix(0.0, 2, 2))
  }

  return(listData)
}

##' @noRd
fitDBVSBVinter_prepIn_setMap <- function(listData, listParams, inMap = NULL) {
  defltMap <- list()
  if (length(listData$y_SC_f) == 1 & length(listData$y_SC_t) == 1) { # only IC
    defltMap <- list(
      beta_SC_f = rep(NA, length(listParams$beta_SC_f)),
      SIGV_f = rep(NA, length(listParams$SIGV_f)),
      log_sd_e_SC_f = NA,
      log_sd_SIGV_f = NA,
      beta_SC_t = rep(NA, length(listParams$beta_SC_t)),
      log_sd_e_SC_t = NA
    )
  } else if (length(listData$y_SC_f) > 1 & length(listData$y_SC_t) == 1) { # IC and SC_f
    defltMap <- list(
      beta_SC_t       = rep(NA, length(listParams$beta_SC_t)),
      log_sd_e_SC_t   = NA
    )
  } else if (length(listData$y_SC_f) == 1 & length(listData$y_SC_t) > 1) { # IC and SC_t
    defltMap <- list(
      beta_SC_f = rep(NA, length(listParams$beta_SC_f)),
      SIGV_f = rep(NA, length(listParams$SIGV_f)),
      log_sd_e_SC_f = NA,
      log_sd_SIGV_f = NA
    )
  }
  if (nrow(listData$Z_DxS_f) <= 2) { # DBVxSBV are ignored
    defltMap$DBVxSBV <- matrix(NA, ncol(listData$Z_DxS_f), 2)
    defltMap$log_sd_DxS <- rep(NA, 2)
    defltMap$unconstr_cor_DxS <- NA
  }
  defltMap <- lapply(defltMap, as.factor)

  if (is.null(inMap)) {
    inMap <- defltMap
  } else {
    stopifnot(!is.null(names(inMap)))
    inMap[intersect(names(inMap), names(defltMap))] <- NULL
    inMap <- c(defltMap, inMap)
  }

  return(inMap)
}

##' @noRd
fitDBVSBVinter_prepIn <- function(listY, listX, listZ, listVCov, REML = TRUE,
                                  lOptions = NULL) {
  out <- list()

  listData <- fitDBVSBVinter_prepIn_listData(listY, listX, listZ, listVCov)
  out$listData <- listData

  ## TODO: get reasonable starting values for genetic (co)variances (fit first with lme4?)
  listParams <- list(
    B_IC = matrix(0.0, ncol(listData$X_IC), 2),
    BV_f = matrix(0.0, ncol(listData$Z_DS_f), 2),
    log_sd_BV_f = rep(log(1), 2),
    unconstr_cor_DS_f = 0.0,
    DBVxSBV = matrix(0.0, ncol(listData$Z_DxS_f), 2),
    log_sd_DxS = rep(log(1), 2),
    unconstr_cor_DxS = 0.0,
    log_sd_E_IC = rep(log(1), 2),
    unconstr_cor_E_IC = 0.0
  )
  listParams$beta_SC_f <- rep(0.0, ncol(listData$X_SC_f))
  listParams$SIGV_f <- rep(0.0, ncol(listData$Z_D_f))
  listParams$log_sd_SIGV_f <- log(1)
  listParams$log_sd_e_SC_f <- log(1)
  listParams$beta_SC_t <- rep(0.0, ncol(listData$X_SC_t))
  listParams$log_sd_e_SC_t <- log(1)
  if (!is.null(lOptions)) {
    if ("listParams" %in% names(lOptions)) {
      listParams <- lOptions$listParams
    }
  }
  out$listParams <- listParams

  vecRnd <- c("BV_f")
  if (nrow(listData$Z_DxS_f) > 2) {
    vecRnd <- c(vecRnd, "DBVxSBV")
  }
  if ("y_SC_f" %in% names(listY)) {
    vecRnd <- c(vecRnd, "SIGV_f")
  }
  if (REML) {
    vecRnd <- c(vecRnd, "B_IC")
    if ("y_SC_f" %in% names(listY)) {
      vecRnd <- c(vecRnd, "beta_SC_f")
    }
    if ("y_SC_t" %in% names(listY)) {
      vecRnd <- c(vecRnd, "beta_SC_t")
    }
  }
  out$vecRnd <- vecRnd

  out$map <- fitDBVSBVinter_prepIn_setMap(
    listData, listParams,
    lOptions[["map"]]
  )

  return(out)
}

##' @noRd
fitDBVSBVinter_prepOut <- function(outTmb, listData, sep = "_") {
  idx <- which(rownames(outTmb$sry_sdr) == "B_IC")
  idx <- idx[1:(ncol(listData$X_IC) * ncol(listData$Y_IC))]
  outTmb$B_IC <- outTmb$sry_sdr[idx, ]
  rownames(outTmb$B_IC) <- paste(
    rownames(outTmb$B_IC),
    rep(colnames(listData$X_IC),
      times = ncol(listData$Y_IC)
    ),
    "on species",
    rep(1:ncol(listData$Y_IC),
      each = ncol(listData$X_IC)
    )
  )

  idx <- which(rownames(outTmb$sry_sdr) %in% c("vars_BV_f", "cor_BV_f"))
  stopifnot(length(idx) == 3)
  outTmb$vcov_BV_f <- outTmb$sry_sdr[idx, ]
  rownames(outTmb$vcov_BV_f)[1:2] <- c("var_DBV_f", "var_SBV_IC_f")

  idx <- which(rownames(outTmb$sry_sdr) == "vars_E_IC")
  stopifnot(length(idx) == ncol(listData$Y_IC))
  outTmb$vcov_E_IC <- outTmb$sry_sdr[idx, ]
  rownames(outTmb$vcov_E_IC) <- paste(
    rownames(outTmb$vcov_E_IC),
    "of species",
    1:ncol(listData$Y_IC)
  )
  idx <- which(rownames(outTmb$sry_sdr) == "Cor_E_IC")
  stopifnot(length(idx) == ncol(listData$Y_IC) * ncol(listData$Y_IC))
  idx_offdiag <- which(upper.tri(matrix(NA, ncol(listData$Y_IC), ncol(listData$Y_IC))))
  idx <- idx[idx_offdiag]
  outTmb$vcov_E_IC <- rbind(
    outTmb$vcov_E_IC,
    outTmb$sry_sdr[idx, ]
  )
  idx <- (ncol(listData$Y_IC) + 1):(nrow(outTmb$vcov_E_IC))
  rownames(outTmb$vcov_E_IC)[idx] <- "cor_E_IC"
  if (ncol(listData$Y_IC) > 2) { # TODO: add "of species i-j"
    rownames(outTmb$vcov_E_IC)[idx] <- paste(
      rownames(outTmb$vcov_E_IC)[idx],
      1:length(idx)
    )
  }

  idx <- which(rownames(outTmb$sry_sdr) == "BV_f")
  stopifnot(length(idx) == nrow(listData$K) * 2)
  idx_DBV <- 1:nrow(listData$K)
  outTmb$DBV_f <- outTmb$sry_sdr[idx[idx_DBV], ]
  rownames(outTmb$DBV_f) <- colnames(listData$Z_DS_f)[idx_DBV]
  idx_SBV <- (nrow(listData$K) + 1):(2 * nrow(listData$K))
  outTmb$SBV_IC_f <- outTmb$sry_sdr[idx[idx_SBV], ]
  rownames(outTmb$SBV_IC_f) <- colnames(listData$Z_DS_f)[idx_DBV]

  idx <- which(rownames(outTmb$sry_sdr) == "BV_IC_f")
  stopifnot(length(idx) == nrow(listData$K))
  outTmb$BV_IC_f <- outTmb$sry_sdr[idx, ]
  rownames(outTmb$BV_IC_f) <- rownames(outTmb$DBV_f)

  if (nrow(listData$Z_DxS_f) <= 2) {
    outTmb$report$DBVxSBV <- NULL
    outTmb$report$log_sd_DxS <- NULL
    outTmb$report$unconstr_cor_DxS <- NULL
  }
  if (length(listData$y_SC_f) == 1 & length(listData$y_SC_t) == 1) { # only IC
    outTmb$report$beta_SC_f <- NULL
    outTmb$report$var_SIGV_f <- NULL
    outTmb$report$var_err_SC_f <- NULL
    outTmb$report$beta_SC_t <- NULL
    outTmb$report$var_err_SC_t <- NULL
    outTmb$report$BV_SC_f <- NULL
  } else if (length(listData$y_SC_f) > 1 & length(listData$y_SC_t) == 1) { # IC and SC_f
    outTmb$report$beta_SC_t <- NULL
    outTmb$report$var_err_SC_t <- NULL
  } else if (length(listData$y_SC_f) == 1 & length(listData$y_SC_t) > 1) { # IC and SC_t
    outTmb$report$beta_SC_f <- NULL
    outTmb$report$var_SIGV_f <- NULL
    outTmb$report$var_err_SC_f <- NULL
    outTmb$report$BV_SC_f <- NULL
  }

  idx <- which(rownames(outTmb$sry_sdr) == "beta_SC_f")
  if (length(idx) > 0) {
    stopifnot(length(idx) == ncol(listData$X_SC_f))
    outTmb$beta_SC_f <- outTmb$sry_sdr[idx, ]
    rownames(outTmb$beta_SC_f) <- paste(
      rownames(outTmb$beta_SC_f),
      colnames(listData$X_SC_f)
    )
  }
  idx <- which(rownames(outTmb$sry_sdr) == "var_SIGV_f")
  if (length(idx) > 0) {
    stopifnot(length(idx) == 1)
    outTmb$var_SIGV_f <- outTmb$sry_sdr[idx, , drop = FALSE]
  }
  idx <- which(rownames(outTmb$sry_sdr) == "SIGV_f")
  if (length(idx) > 0) {
    stopifnot(length(idx) == ncol(listData$Z_D_f))
    outTmb$SIGV_f <- outTmb$sry_sdr[idx, ]
    rownames(outTmb$SIGV_f) <- colnames(listData$Z_D_f)
  }
  idx <- which(rownames(outTmb$sry_sdr) == "BV_SC_f")
  if (length(idx) > 0) {
    stopifnot(length(idx) == ncol(listData$Z_D_f))
    outTmb$BV_SC_f <- outTmb$sry_sdr[idx, ]
    rownames(outTmb$BV_SC_f) <- colnames(listData$Z_D_f)
  }
  idx <- which(rownames(outTmb$sry_sdr) == "var_err_SC_f")
  if (length(idx) > 0) {
    stopifnot(length(idx) == 1)
    outTmb$var_err_SC_f <- outTmb$sry_sdr[idx, , drop = FALSE]
  }

  idx <- which(rownames(outTmb$sry_sdr) == "beta_SC_t")
  if (length(idx) > 0) {
    stopifnot(length(idx) == ncol(listData$X_SC_t))
    outTmb$beta_SC_t <- outTmb$sry_sdr[idx, ]
    rownames(outTmb$beta_SC_t) <- paste(
      rownames(outTmb$beta_SC_t),
      colnames(listData$X_SC_t)
    )
  }
  idx <- which(rownames(outTmb$sry_sdr) == "var_err_SC_t")
  if (length(idx) > 0) {
    stopifnot(length(idx) == 1)
    outTmb$var_err_SC_t <- outTmb$sry_sdr[idx, , drop = FALSE]
  }

  return(outTmb)
}

##' @noRd
fitDBVSBVinter_AIC <- function(outTmb, listData) {
  nbPars <- prod(dim(outTmb$report$B_IC)) +
    (nrow(outTmb$report$Cor_BV_f) * (nrow(outTmb$report$Cor_BV_f) + 1)) / 2 +
    (nrow(outTmb$report$Cor_E_IC) * (nrow(outTmb$report$Cor_E_IC) + 1)) / 2
  if (nrow(listData$Z_DxS_f) > 2) {
    nbPars <- nbPars +
      (nrow(outTmb$report$Cor_DxS) * (nrow(outTmb$report$Cor_DxS) + 1)) / 2
  }
  if (length(listData$y_SC_f) > 1) {
    nbPars <- nbPars +
      length(outTmb$report$beta_SC_f) +
      1 + # var(SIGV_f)
      1 # var(err_SC_f)
  }
  if (length(listData$y_SC_t) > 1) {
    nbPars <- nbPars +
      length(outTmb$report$beta_SC_t) +
      1 # var(err_SC_t)
  }
  infoCriterion(
    k = nbPars,
    lnLmax = -outTmb$fit$objective, # "-" because the *neg*loglik was maximized
    type = "AIC"
  )
}

##' Fit DBV-SBV models for interspecific mixtures
##'
##' Fits DBV-SBV models for interspecific mixtures.
##' @param listY named list of response variables (a two-column matrix Y_IC, with the first column for the focal species and the second column for the tester species in the case of a tester-based design, and, optionally, a vector y_SC_f and a vector y_SC_t)
##' @param listX named list of design matrices for the fixed-effect explanatory factors (a matrix X_IC and, optionally, a matrix X_SC_f and a matrix X_SC_t)
##' @param listZ named list of design matrices of random-effect explanatory factors (a matrix Z_DS_f and, optionally, a matrix Z_D_f)
##' @param listVCov named list of square, symmetric matrices used, after re-scaling, as variance-covariance matrices for the random-effect factors; named "K" for the BVs (DBVs and SBVs) and "Kmixpair" for the DBVxSBVs; these matrices must have dimnames (both rows and columns) coherent with the column names of the design matrices in \code{listZ}
##' @param REML logical
##' @param lOptions named list of options (for experts)
##' @param verbose verbosity level
##' @return list
##' @seealso \code{\link{simulDBVSBVinter}}
##' @author Jemay Salomon, Timothee Flutre
##' @examples
##' ## simulate a data set with both sole crops and intercrops:
##' GRMs <- list("S1" = diag(100),
##'              "S2" = diag(2))
##' dimnames(GRMs$S1) <- list(paste0("gS1_", 1:100), paste0("gS1_", 1:100))
##' dimnames(GRMs$S2) <- list(paste0("gS2_", 1:2), paste0("gS2_", 1:2))
##' out <- simulDBVSBVinter(GRMs)
##' names(out)
##' datW <- out$datW
##' str(datW)
##'
##' ## fit the model using intercrops only:
##' idxIC <- which(!is.na(datW$geno_S1) & !is.na(datW$geno_S2))
##' datW_IC <- droplevels(datW[idxIC, ])
##' listY <- list(Y_IC = datW_IC[, c("yield_S1", "yield_S2")])
##' listX <- list(X_IC = model.matrix(~ 1 + block + geno_S2, datW_IC,
##'               contrasts.arg = list("block" = "contr.sum",
##'                                    "geno_S2" = "contr.sum")))
##' listZ <- list(Z_DS_f = model.matrix(~ 0 + geno_S1, datW_IC))
##' colnames(listZ$Z_DS_f) <- gsub("^geno_S1", "", colnames(listZ$Z_DS_f))
##' listVCov <- list(K = GRMs$S1[levels(datW_IC$geno_S1), levels(datW_IC$geno_S1)])
##' fitTmb <- fitDBVSBVinter(listY, listX, listZ, listVCov,
##'                          lOptions = list(iter.max = 20), REML = TRUE)
##' names(fitTmb)
##' fitTmb$sdr
##'
##' ## see the third vignette for more details
##' @export
fitDBVSBVinter <- function(listY, listX, listZ, listVCov, REML = TRUE, lOptions = NULL,
                           verbose = FALSE) {
  ## checks
  stopifnot(
    is.list(listY),
    length(listY) <= 3,
    !is.null(names(listY)),
    all(names(listY) %in% c("Y_IC", "y_SC_f", "y_SC_t")),
    "Y_IC" %in% names(listY),
    ncol(listY$Y_IC) == 2,
    is.list(listX),
    length(listX) <= 3,
    !is.null(names(listX)),
    all(names(listX) %in% c("X_IC", "X_SC_f", "X_SC_t")),
    "X_IC" %in% names(listX),
    is.list(listZ),
    all(sapply(listZ, is.matrix)),
    !is.null(names(listZ)),
    all(names(listZ) %in% c("Z_DS_f", "Z_D_f", "Z_DxS_f")),
    "Z_DS_f" %in% names(listZ),
    is.list(listVCov),
    all(sapply(listVCov, is.matrix)),
    !is.null(names(listVCov)),
    "K" %in% names(listVCov),
    nrow(listVCov$K) == ncol(listVCov$K),
    !is.null(dimnames(listVCov$K)),
    all(rownames(listVCov$K) == colnames(listVCov$K)),
    all(colnames(listZ$Z_DS_f) %in% rownames(listVCov$K)),
    all(colnames(listZ$Z_DS_f) == rep(rownames(listVCov$K), times = 2)),
    is.logical(REML)
  )
  if ("Kmixpair" %in% names(listVCov)) {
    stopifnot(
      nrow(listVCov$Kmixpair) == ncol(listVCov$Kmixpair),
      !is.null(dimnames(listVCov$Kmixpair)),
      all(rownames(listVCov$Kmixpair) == colnames(listVCov$Kmixpair)),
      "Z_DxS_f" %in% names(listZ),
      all(colnames(listZ$Z_DxS_f) %in% rownames(listVCov$Kmixpair))
    )
  }
  if (length(listY) == 2) {
    stopifnot(
      "y_SC_f" %in% names(listY),
      length(listX) == 2,
      "X_SC_f" %in% names(listX),
      length(listZ) == 2,
      "Z_D_f" %in% names(listZ)
    )
  }
  if (length(listY) == 3) {
    stopifnot(
      "y_SC_t" %in% names(listY),
      length(listX) == 3,
      "X_SC_t" %in% names(listX)
    )
  }
  if (!is.null(lOptions)) {
    stopifnot(
      is.list(lOptions),
      !is.null(names(lOptions)),
      all(names(lOptions) %in% c(
        "listParams", "map", "makeAD", "optim",
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
  inputs4TMB <- fitDBVSBVinter_prepIn(listY, listX, listZ, listVCov, REML, lOptions)
  out$inputs4TMB <- inputs4TMB[c("listParams", "map", "vecRnd")]
  if (verbose) {
    str(inputs4TMB)
  }

  if (lOptions$makeAD) {
    if (verbose) {
      print("automatic differentiation")
    }
    f <- MakeADFun(
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
  }

  if (lOptions$makeAD & lOptions$optim) {
    if (verbose) {
      print("optimization")
    }
    capture <- capture.output(
      fit <- nlminb(
        start = f$par, objective = f$fn, gradient = f$gr, hessian = NULL,
        control = list(
          "trace" = 1,
          "iter.max" = lOptions$iter.max
        )
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
    out$sry_sdr <- summary(out$sdr, select = "report", p.value = TRUE)

    if (verbose) {
      print("output preparation")
    }
    out <- fitDBVSBVinter_prepOut(out, inputs4TMB$listData, "_")
    if (FALSE) {
      out$AIC <- fitDBVSBVinter_AIC(out, inputs4TMB$listData)
    }
  }

  return(out)
}
