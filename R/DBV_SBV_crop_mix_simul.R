##' Create the data structure including all sole crops and all possible intercrops as a `data.frame` for one block
##' @noRd
makeCompleteDesignOneBlockInterMixDesign <- function(levGenos, sep = "+") {
  levSpecies <- names(levGenos)
  S1 <- names(levGenos)[1]
  S2 <- names(levGenos)[2]
  nbGenos <- sapply(levGenos, length)
  dat_bl1 <- data.frame(
    standID = NA,
    geno_S1 = c(
      levGenos[[S1]], # each geno of species 1 as sole crop
      rep(NA, nbGenos[S2]), # each geno of species 2 as sole crop
      rep(levGenos[[S1]], each = nbGenos[S2])
    ), # all intercrops
    geno_S2 = c(
      rep(NA, nbGenos[S1]),
      levGenos[[S2]],
      rep(levGenos[[S2]], nbGenos[S1])
    ),
    type = c(
      rep("SC", nbGenos[S1] + nbGenos[S2]),
      rep("IC", nbGenos[S1] * nbGenos[S2])
    ),
    type2 = c(
      rep("sole_S1", nbGenos[S1]),
      rep("sole_S2", nbGenos[S2]),
      rep("mix", nbGenos[S1] * nbGenos[S2])
    ),
    block = NA,
    x = NA,
    y = NA,
    plot = NA,
    true_gen_yield_S1 = NA, # only genetic factors
    true_gen_yield_S2 = NA, # idem
    true_yield_S1 = NA, # true_gen_yield + mu + block factors
    true_yield_S2 = NA, # idem
    true_fix_yield_S1 = NA, # only fixed effects
    true_fix_yield_S2 = NA, # only fixed effects
    true_rnd_yield_S1 = NA, # only random effects
    true_rnd_yield_S2 = NA, # only random effects
    yield_S1 = NA, # true_yield + error
    yield_S2 = NA
  )

  dat_bl1$standID <- paste0(dat_bl1$geno_S1, sep, dat_bl1$geno_S2)
  idx <- which(is.na(dat_bl1$geno_S2))
  dat_bl1$standID[idx] <- dat_bl1$geno_S1[idx]
  idx <- which(is.na(dat_bl1$geno_S1))
  dat_bl1$standID[idx] <- dat_bl1$geno_S2[idx]
  return(dat_bl1)
}

##' Create the data structure including all sole crops and all possible intercrops as a `data.frame` for all blocks
##' @noRd
makeCompleteDesignInterMixDesign <- function(dat_bl1, levBlocks) {
  datW <- dat_bl1
  nbBlocks <- length(levBlocks)
  if (nbBlocks > 1) {
    for (k in 2:nbBlocks) {
      datW <- rbind(datW, dat_bl1)
    }
  }
  datW$standID <- factor(datW$standID)
  datW$geno_S1 <- factor(datW$geno_S1)
  datW$geno_S2 <- factor(datW$geno_S2)
  datW$type <- factor(datW$type, levels = c("SC", "IC"))
  datW$type2 <- factor(datW$type2, levels = c("sole_S1", "sole_S2", "mix"))
  datW$block <- factor(rep(levBlocks, each = nrow(dat_bl1)))
  return(datW)
}

##' Make a balanced, incomplete design by sampling a subset of sole crops and intercrops
##' @noRd
allocGenosMixsInterMixDesign <- function(datW, nbBlocks, levGenos) {
  nbGenos <- sapply(levGenos, length)

  ## all genos of species 1 obs. as sole crop:
  ## 50% in block A and 50% in block B
  soleGenosA <- sample(levGenos[["S1"]], nbGenos["S1"] / nbBlocks)
  soleGenosB <- levGenos[["S1"]][!levGenos[["S1"]] %in% soleGenosA]
  stopifnot(
    anyDuplicated(c(soleGenosA, soleGenosB)) == 0,
    length(c(soleGenosA, soleGenosB)) == nbGenos["S1"]
  )

  ## 50% obs with tester 1 and 50% obs with tester 2:
  ## 50% of those in sole crop in A will be in intercrop with the 1st tester in A
  ## 50% of those in sole crop in A will be in intercrop with the 2nd tester in B
  ## 50% of those in sole crop in B will be in intercrop with the 1st tester in B
  ## 50% of those in sole crop in B will be in intercrop with the 2nd tester in A
  interGenos <- list(half1A = sample(soleGenosA, length(soleGenosA) / 2))
  interGenos$half1B <- sample(
    soleGenosA[!soleGenosA %in% interGenos$half1A],
    length(soleGenosA) / 2
  )
  interGenos$half2A <- sample(soleGenosB, length(soleGenosB) / 2)
  interGenos$half2B <- sample(
    soleGenosB[!soleGenosB %in% interGenos$half2A],
    length(soleGenosB) / 2
  )
  stopifnot(
    all(sort(do.call(c, interGenos[grep("1", names(interGenos))])) == sort(soleGenosA)),
    all(sort(do.call(c, interGenos[grep("2", names(interGenos))])) == sort(soleGenosB))
  )
  stopifnot(
    names(table(table(do.call(c, interGenos)))) == 1,
    length(do.call(c, interGenos)) == nbGenos["S1"]
  )

  ## subset the complete design:
  idxSole <- c(
    which(as.character(datW$geno_S1) %in% soleGenosA &
      is.na(datW$geno_S2) & datW$block == "A"),
    which(as.character(datW$geno_S1) %in% soleGenosB &
      is.na(datW$geno_S2) & datW$block == "B")
  )
  stopifnot(length(idxSole) == nbGenos["S1"])
  idxInter <- c(
    which(datW$geno_S1 %in% interGenos$half1A &
      datW$geno_S2 == levGenos[["S2"]][1] &
      datW$block == "A"),
    which(datW$geno_S1 %in% interGenos$half2B &
      datW$geno_S2 == levGenos[["S2"]][2] &
      datW$block == "B"),
    which(datW$geno_S1 %in% interGenos$half1B &
      datW$geno_S2 == levGenos[["S2"]][1] & datW$block == "B"),
    which(datW$geno_S1 %in% interGenos$half2A &
      datW$geno_S2 == levGenos[["S2"]][2] &
      datW$block == "A")
  )
  stopifnot(length(idxInter) == nbGenos["S1"])

  ## keep the sole crop of species 2:
  idx2 <- which(datW$type2 == "sole_S2")

  datW <- droplevels(datW[c(idxSole, idxInter, idx2), ])
  rownames(datW) <- NULL
  return(datW)
}

##' @noRd
plotAllocSchemeInterMixDesign <- function(datW, title = "Sparse trial design for both sole crops and intercrops") {
  idx <- which(datW$type2 == "sole_S2")
  tmp <- droplevels(datW[-idx, c("geno_S1", "geno_S2", "type2", "block")])
  ## envs for wheat:
  ##  sole crops: SC_A, SC_B
  ##  intercrops: IC_P1_A, IC_P2_A, IC_P2_A, IC_P2_B
  tmp$env <- ""
  tmp$env[grep("^sole", tmp$type2)] <- "sole crop"
  tmp$env[which(tmp$type2 == "mix")] <- "intercrop"
  for (i in 1:nlevels(tmp$geno_S2)) {
    idx <- which(tmp$type2 == "mix" & tmp$geno_S2 == levels(tmp$geno_S2)[i])
    tmp$env[idx] <- paste0(tmp$env[idx], "\nwith pea", i)
  }
  tmp$env <- paste0(tmp$env, "\nblock ", tmp$block)
  tmp$env <- factor(tmp$env,
    levels = rev(c(
      paste0("sole crop\nblock ", levels(tmp$block)),
      paste0(
        rep(paste0("intercrop\nwith pea", 1:nlevels(tmp$geno_S2), "\nblock "),
          each = nlevels(tmp$block)
        ),
        levels(tmp$block)
      )
    ))
  )

  p <- ggplot(tmp) +
    aes(x = .data[["geno_S1"]], y = .data[["env"]], fill = .data[["env"]]) +
    geom_tile(show.legend = FALSE) +
    labs(
      title = title,
      subtitle = paste0(
        nrow(tmp), " micro-plots",
        "; ", nlevels(tmp$geno_S1), " focal genotypes (species S1)",
        "; ", nlevels(tmp$geno_S2), " tester genotypes (species S2)",
        "; ", nlevels(tmp$block), " randomized incomplete blocks"
      ),
      x = "focal genotypes",
      y = ""
    ) + # "wheat environments") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  return(p)
}

##' @noRd
setMplotCoordsInterMixDesign <- function(datW, nbYs = 10) {
  datW$plot <- factor(paste0(
    as.character(datW$geno_S1),
    "_", as.character(datW$geno_S2),
    "_", as.character(datW$block)
  ))
  stopifnot(nlevels(datW$plot) == nrow(datW))
  nbXs <- ceiling(nlevels(datW$plot) / nbYs)
  stopifnot(nbXs * nbYs >= nlevels(datW$plot))
  seq_y <- seq(1, nbYs)
  blockWidth <- ceiling(nbXs / nlevels(datW$block))
  for (i in 1:nlevels(datW$block)) {
    block <- levels(datW$block)[i]
    plots <- unique(as.character(datW$plot)[datW$block == block])
    first_x <- (i - 1) * blockWidth + 1
    last_x <- first_x + blockWidth - 1
    seq_x <- seq(first_x, last_x)
    coords <- data.frame(
      x = rep(seq_x, each = nbYs),
      y = seq_y
    )
    stopifnot(nrow(coords) >= length(plots))
    idxPlotCoords <- sample(1:nrow(coords), length(plots))
    for (j in seq_along(plots)) {
      idx <- which(datW$plot == plots[j])
      datW[idx, c("x", "y")] <- coords[idxPlotCoords[j], ]
    }
  }
  stopifnot(
    all(!is.na(datW$x)),
    all(!is.na(datW$y))
  )
  datW$plot <- factor(paste0(datW$x, as.character(datW$block), datW$y))
  return(datW)
}

##' @noRd
drawParamsInterMixDesign <- function(datW, GRMs, sep,
                                     parsGenetics = list(
                                       mu = list(
                                         "S1" = c("SC" = 65, "IC" = 32),
                                         "S2" = c("SC" = 30, "IC" = 27)
                                       ),
                                       sdBlocks = 2,
                                       CV_g = c("S1" = 0.08, "S2" = 0.08),
                                       prop_var_SBV = c("S1" = 0.2, "S2" = 0.2),
                                       prop_var_SIGV = c("S1" = 0.2, "S2" = 0),
                                       cor_DBV_SBV = c("S1" = -0.9, "S2" = -0.9),
                                       prop_var_DBVxSBV = c("S1" = 0, "S2" = 0),
                                       use_GRM_for_DBVxSBV = TRUE,
                                       cor_DxS = 0,
                                       H2_SC = c("S1" = 0.7, "S2" = 0.7),
                                       T2 = c("S1" = 0.7, "S2" = 0.7),
                                       cor_err_IC = -0.2
                                     ),
                                     parsDesign = list(
                                       nbBlocks = 2,
                                       nbPlotsPerBl = 500,
                                       nbYs = 10,
                                       testersS2 = TRUE
                                     )) {
  mu <- NA
  sdBlocks <- NA
  CV_g <- NA
  prop_var_SBV <- NA
  prop_var_SIGV <- NA
  cor_DBV_SBV <- NA
  prop_var_DBVxSBV <- NA
  use_GRM_for_DBVxSBV <- NA
  cor_DxS <- NA
  H2_SC <- NA
  T2 <- NA
  cor_err_IC <- NA # avoid NOTE in R CMD check
  e <- list2env(parsGenetics, envir = environment())

  levSpecies <- names(mu)
  levBlocks <- levels(datW$block)
  nbBlocks <- length(levBlocks)
  levGenos <- lapply(GRMs, rownames)
  nbGenos <- sapply(levGenos, length)

  truth <- parsGenetics

  ## effects of the "block" factor
  blEffs <- list()
  for (species in levSpecies) {
    tmp <- rnorm(n = nbBlocks, mean = 0, sd = sdBlocks)
    tmp <- as.vector(scale(tmp, center = TRUE, scale = FALSE))
    blEffs[[species]] <- setNames(tmp, levBlocks)
  }
  truth[["blEffs"]] <- blEffs

  ## true means of the "block" factor
  blTrueMeans <- list()
  for (species in levSpecies) {
    tmp <- matrix(
      nrow = nbBlocks, ncol = 2,
      dimnames = list(levBlocks, c("SC", "IC"))
    )
    for (type in c("SC", "IC")) {
      for (block in levBlocks) {
        tmp[block, type] <- mu[[species]][type] + blEffs[[species]][block]
      }
    }
    blTrueMeans[[species]] <- tmp
  }
  truth[["blTrueMeans"]] <- blTrueMeans

  ## true contrasts of the "block" factor
  blTrueContrs <- list()
  if (nbBlocks > 1) {
    Bstar.cs <- contr.sum(levels(datW$block))
    B.cs <- cbind(intercept = 1, Bstar.cs)
    C.cs <- solve(B.cs)
    blTrueContrs$S1 <- C.cs %*% blTrueMeans$S1
    blTrueContrs$S2 <- C.cs %*% blTrueMeans$S2
  }
  truth[["blTrueContrs"]] <- blTrueContrs

  var_DBV <- c(
    "S1" = as.vector((CV_g["S1"] * mu$S1["SC"])^2),
    "S2" = as.vector((CV_g["S2"] * mu$S2["SC"])^2)
  )
  truth[["var_DBV"]] <- var_DBV

  var_SBV <- prop_var_SBV * var_DBV
  truth[["var_SBV"]] <- var_SBV

  var_SIGV <- prop_var_SIGV * var_DBV
  truth[["var_SIGV"]] <- var_SIGV

  cov_DBV_SBV <- cor_DBV_SBV * sqrt(var_DBV) * sqrt(var_SBV)
  truth[["cov_DBV_SBV"]] <- cov_DBV_SBV

  Sigma_DBV_SBV <- lapply(levSpecies, function(species) {
    matrix(
      data = c(
        var_DBV[species], cov_DBV_SBV[species],
        cov_DBV_SBV[species], var_SBV[species]
      ),
      nrow = 2, ncol = 2, byrow = TRUE,
      dimnames = list(
        c("DBV", "SBV"),
        c("DBV", "SBV")
      )
    )
  })
  names(Sigma_DBV_SBV) <- levSpecies
  truth[["Sigma_DBV_SBV"]] <- Sigma_DBV_SBV

  var_DBVxSBV <- prop_var_DBVxSBV * var_DBV
  truth[["var_DBVxSBV"]] <- var_DBVxSBV

  cov_DBV_SIGV <- 0 # TODO: allow this to be non-zero
  var_cBV <- var_DBV + var_SIGV + cov_DBV_SIGV
  truth[["var_cBV"]] <- var_cBV

  var_TBV <- var_DBV + 2 * cov_DBV_SBV + var_SBV # full dilution: eq.5 from Bijma (2010)
  ## TODO: if var(DBVxSBV) > 0, need to take it into account in T2?
  truth[["var_TBV"]] <- var_TBV

  BV <- lapply(levSpecies, function(species) {
    M <- matrix(0,
      nrow = nbGenos[species], ncol = 2,
      dimnames = list(levGenos[[species]], colnames(Sigma_DBV_SBV[[species]]))
    )
    rmatnorm(
      n = 1, M = M,
      U = GRMs[[species]][
        levGenos[[species]],
        levGenos[[species]]
      ],
      V = Sigma_DBV_SBV[[species]]
    )[, , 1]
  })
  names(BV) <- levSpecies
  if (parsDesign$testersS2) {
    BV$S2[, "DBV"] <- scale(BV$S2[, "DBV"], center = TRUE, scale = FALSE)
    BV$S2[, "SBV"] <- scale(BV$S2[, "SBV"], center = TRUE, scale = FALSE)
  }
  truth[["BV"]] <- BV
  truth[["BV_IC"]] <- lapply(BV, rowSums)

  BVTrueMeans <- list()
  for (species in levSpecies) {
    otherSpecies <- levSpecies[levSpecies != species]
    tmp <- array(
      dim = c(nlevels(datW[[paste0("geno_", species)]]), 2, 2),
      dimnames = list(
        levels(datW[[paste0("geno_", species)]]),
        c("SC", "IC"),
        c("DBV", "SBV")
      )
    )
    for (geno in dimnames(tmp)[[1]]) {
      for (type in dimnames(tmp)[[2]]) {
        tBV <- "DBV"
        tmp[geno, type, tBV] <- mu[[species]][type] + BV[[species]][geno, tBV]
        if (type == "IC") {
          tBV <- "SBV"
          tmp[geno, type, tBV] <- mu[[otherSpecies]][type] + BV[[species]][geno, tBV]
        }
      }
    }
    BVTrueMeans[[species]] <- tmp
  }
  truth[["BVTrueMeans"]] <- BVTrueMeans

  BVTrueContrs <- list()
  idxIC <- which(!is.na(datW$geno_S1) & !is.na(datW$geno_S2))
  for (species in levSpecies) {
    Bstar.cs <- contr.sum(levels(datW[[paste0("geno_", species)]][idxIC]))
    B.cs <- cbind(intercept = 1, Bstar.cs)
    C.cs <- solve(B.cs)
    BVTrueContrs[[species]] <- list(
      SC_DBV = C.cs %*% BVTrueMeans[[species]][, "SC", "DBV"],
      IC_DBV = C.cs %*% BVTrueMeans[[species]][, "IC", "DBV"],
      IC_SBV = C.cs %*% BVTrueMeans[[species]][, "IC", "SBV"]
    )
  }
  truth[["BVTrueContrs"]] <- BVTrueContrs

  SIGVs <- lapply(levSpecies, function(species) {
    tmp <- mvrnorm(
      n = 1, mu = rep(0, nbGenos[species]),
      Sigma = var_SIGV[species] * GRMs[[species]]
    )
    names(tmp) <- levGenos[[species]]
    tmp
  })
  names(SIGVs) <- levSpecies
  if (parsDesign$testersS2) {
    tmp <- as.vector(scale(SIGVs$S2, center = TRUE, scale = FALSE))
    SIGVs$S2 <- setNames(tmp, names(SIGVs$S2))
  }
  truth[["SIGVs"]] <- SIGVs

  BV_SC <- lapply(levSpecies, function(species) {
    stopifnot(all(rownames(BV[[species]]) == names(SIGVs[[species]])))
    BV[[species]][, "DBV"] + SIGVs[[species]]
  })
  names(BV_SC) <- levSpecies
  truth[["BV_SC"]] <- BV_SC

  ## Create a matrix of dimensions (nbGenosS1 x nbGenosxS2) x 2
  ## * the rows represent the intercrops
  ## * the columns represent the species
  ## It extends eq.29 in suppl section 1.2 of ATOMM (Wang et al, 2018)
  vcovMatToUse <- diag(prod(nbGenos))
  if (use_GRM_for_DBVxSBV) {
    vcovMatToUse <- GRMs[[2]] %x% GRMs[[1]]
    rownames(vcovMatToUse) <- paste0(
      rep(levGenos[[2]],
        each = nbGenos[1]
      ),
      sep,
      levGenos[[1]]
    )
    colnames(vcovMatToUse) <- rownames(vcovMatToUse)
  }
  M <- matrix(0, nrow = prod(nbGenos), ncol = 2)
  rownames(M) <- rownames(vcovMatToUse)
  colnames(M) <- levSpecies
  Sigma_DxS <- matrix(NA, nrow = 2, ncol = 2, dimnames = list(levSpecies, levSpecies))
  Sigma_DxS[1, 1] <- prop_var_DBVxSBV[1] * var_DBV[1]
  Sigma_DxS[2, 2] <- prop_var_DBVxSBV[2] * var_DBV[2]
  Sigma_DxS[1, 2] <- Sigma_DxS[2, 1] <- cor_DxS * sqrt(prod(diag(Sigma_DxS)))
  DBVxSBV <- M
  if (all(diag(Sigma_DxS) != 0)) {
    DBVxSBV <- rmatnorm(n = 1, M = M, U = vcovMatToUse, V = Sigma_DxS)[, , 1]
  }
  ## TODO: use scale() for S2?
  truth[["DBVxSBV"]] <- DBVxSBV

  var_err_SC <- ((1 - H2_SC) / H2_SC) * var_cBV # ignore nb of reps because incomplete design
  truth[["var_err_SC"]] <- var_err_SC
  var_err_IC <- ((1 - T2) / T2) * var_TBV
  truth[["var_err_IC"]] <- var_err_IC
  cov_err_IC <- cor_err_IC * sqrt(var_err_IC["S1"]) * sqrt(var_err_IC["S2"])
  truth[["cov_err_IC"]] <- cov_err_IC
  Sigma_err_IC <- matrix(
    data = c(
      var_err_IC["S1"], cov_err_IC,
      cov_err_IC, var_err_IC["S2"]
    ),
    nrow = 2, ncol = 2,
    dimnames = list(levSpecies, levSpecies)
  )
  truth[["Sigma_err_IC"]] <- Sigma_err_IC

  allErrors <- list()
  for (idx in 1:nrow(datW)) {
    g1 <- as.character(datW$geno_S1[idx])
    g2 <- as.character(datW$geno_S2[idx])
    if (!is.na(g1) & is.na(g2)) { # sole crop of species 1
      species <- levSpecies[1]
      err <- rnorm(n = 1, mean = 0, sd = sqrt(var_err_SC[species]))
    } else if (is.na(g1) & !is.na(g2)) { # sole crop of species 2
      species <- levSpecies[2]
      err <- rnorm(n = 1, mean = 0, sd = sqrt(var_err_SC[species]))
    } else if (!is.na(g1) & !is.na(g2)) { # intercrop
      err <- mvrnorm(n = 1, mu = rep(0, 2), Sigma = Sigma_err_IC)
    }
    allErrors[[idx]] <- err
  }
  truth[["allErrors"]] <- allErrors

  ## fixed-effect matrix and vectors used below to compute yields with
  ## matrix-based formulas:
  B_IC <- NULL
  beta_SC_f <- NULL
  beta_SC_t <- NULL
  if (nbBlocks > 1) {
    if (parsDesign$testersS2) {
      B_IC <- matrix(
        nrow = 1 + nlevels(datW$block) - 1 +
          nlevels(datW$geno_S2) - 1, # should correspond to ncol(X_IC)
        ncol = 2,
        dimnames = list(
          c(
            "(Intercept)",
            paste0("block", 1:(nlevels(datW$block) - 1)),
            paste0("geno_S2", 1:(nlevels(datW$geno_S2) - 1))
          ), # should correspond to colnames(X_IC)
          levSpecies
        )
      )
      for (species in levSpecies) {
        B_IC[1:nrow(blTrueContrs[[species]]), species] <- blTrueContrs[[species]][, "IC"]
      }
      B_IC[(1 + nlevels(datW$block) - 1 + 1):nrow(B_IC), "S1"] <- BVTrueContrs$S2$IC_SBV[-1]
      B_IC[(1 + nlevels(datW$block) - 1 + 1):nrow(B_IC), "S2"] <- BVTrueContrs$S2$IC_DBV[-1]
      beta_SC_f <- blTrueContrs[["S1"]][, "SC"]
      beta_SC_t <- c(
        blTrueContrs[["S2"]][, "SC"],
        BVTrueContrs$S2$SC_DBV[-1, 1]
      )
    } else {
      stop("not yet implemented")
    }
  } else { # nbBlocks == 1
    if (parsDesign$testersS2) {
      B_IC <- matrix(
        nrow = 1 + nlevels(datW$geno_S2) - 1,
        ncol = 2,
        dimnames = list(
          c(
            "(Intercept)",
            paste0("geno_S2", 1:(nlevels(datW$geno_S2) - 1))
          ),
          levSpecies
        )
      )
      for (species in levSpecies) {
        B_IC[1, species] <- mu[[species]]["IC"]
      }
      B_IC[2:nrow(B_IC), "S1"] <- BVTrueContrs$S2$IC_SBV[-1]
      B_IC[2:nrow(B_IC), "S2"] <- BVTrueContrs$S2$IC_DBV[-1]
      beta_SC_f <- mu[["S1"]]["SC"]
      beta_SC_t <- BVTrueContrs$S2$SC_DBV # TODO: + SIGVs[["S2"]]
    } else {
      stop("not yet implemented")
    }
  }
  truth[["B_IC"]] <- B_IC
  truth[["beta_SC_f"]] <- beta_SC_f
  truth[["beta_SC_t"]] <- beta_SC_t

  return(truth)
}

##' @noRd
computeYieldForLoopInterMixDesign <- function(datW, truth, sep) {
  blTrueContrs <- NA
  BV <- NA
  SIGVs <- NA
  mu <- NA
  blEffs <- NA
  allErrors <- NA
  DBVxSBV <- NA # avoid NOTE in R CMD check
  e <- list2env(truth, envir = environment())
  levSpecies <- names(BV)

  deltas_IC_SIGV <- c("S1" = 0, "S2" = 0)
  ## fixed at 0 here
  ## TODO: if != 0 -> c.f. model 3 from Forst et al (2019)

  for (idx in 1:nrow(datW)) {
    true_gen_y1 <- NA
    true_gen_y2 <- NA
    true_y1 <- NA
    true_y2 <- NA
    true_fix_y1 <- NA
    true_fix_y2 <- NA
    true_rnd_y1 <- NA
    true_rnd_y2 <- NA
    y1 <- NA
    y2 <- NA
    g1 <- as.character(datW$geno_S1[idx])
    g2 <- as.character(datW$geno_S2[idx])
    block <- as.character(datW$block[idx])

    if (!is.na(g1) & is.na(g2)) { # sole crop of species 1 (focal)
      species <- levSpecies[1]
      true_gen_y1 <- BV[[species]][g1, "DBV"] + SIGVs[[species]][g1]
      true_y1 <- true_gen_y1 + mu[[species]]["SC"] + blEffs[[species]][block]
      true_fix_y1 <- mu[[species]]["SC"] + blEffs[[species]][block]
      true_rnd_y1 <- true_gen_y1
      y1 <- true_y1 + allErrors[[idx]]
    } else if (is.na(g1) & !is.na(g2)) { # sole crop of species 2 (tester)
      species <- levSpecies[2]
      true_gen_y2 <- BV[[species]][g2, "DBV"] #+ SIGVs[[species]][g2]
      ## SIGVs$S2 being ignored -> see the matrix products below
      true_y2 <- true_gen_y2 + mu[[species]]["SC"] + blEffs[[species]][block]
      true_fix_y2 <- true_y2
      y2 <- true_y2 + allErrors[[idx]]
    } else if (!is.na(g1) & !is.na(g2)) { # intercrop
      S1 <- levSpecies[1]
      S2 <- levSpecies[2]
      true_gen_y1 <- BV[[S1]][g1, "DBV"] + BV[[S2]][g2, "SBV"] +
        DBVxSBV[paste0(g2, sep, g1), S1] +
        deltas_IC_SIGV[[S1]] * SIGVs[[S1]][g1]
      true_y1 <- true_gen_y1 + mu[[S1]]["IC"] + blEffs[[S1]][block]
      true_fix_y1 <- mu[[S1]]["IC"] + blEffs[[S1]][block] + BV[[S2]][g2, "SBV"]
      true_rnd_y1 <- BV[[S1]][g1, "DBV"] +
        DBVxSBV[paste0(g2, sep, g1), S1] +
        deltas_IC_SIGV[[S1]] * SIGVs[[S1]][g1]
      y1 <- true_y1 + allErrors[[idx]][S1]
      true_gen_y2 <- BV[[S2]][g2, "DBV"] + BV[[S1]][g1, "SBV"] +
        DBVxSBV[paste0(g2, sep, g1), S2] +
        deltas_IC_SIGV[[S2]] * SIGVs[[S2]][g2]
      true_y2 <- true_gen_y2 + mu[[S2]]["IC"] + blEffs[[S2]][block]
      true_fix_y2 <- mu[[S2]]["IC"] + blEffs[[S2]][block] + BV[[S2]][g2, "DBV"] +
        deltas_IC_SIGV[[S2]] * SIGVs[[S2]][g2]
      true_rnd_y2 <- BV[[S1]][g1, "SBV"] +
        DBVxSBV[paste0(g2, sep, g1), S2]
      y2 <- true_y2 + allErrors[[idx]][S2]
    }

    datW$true_gen_yield_S1[idx] <- true_gen_y1
    datW$true_gen_yield_S2[idx] <- true_gen_y2
    datW$true_fix_yield_S1[idx] <- true_fix_y1
    datW$true_fix_yield_S2[idx] <- true_fix_y2
    datW$true_rnd_yield_S1[idx] <- true_rnd_y1
    datW$true_rnd_yield_S2[idx] <- true_rnd_y2
    datW$true_yield_S1[idx] <- true_y1
    datW$true_yield_S2[idx] <- true_y2
    datW$yield_S1[idx] <- y1
    datW$yield_S2[idx] <- y2
  }

  return(datW)
}

computeObsMeansContrsInterMixDesign <- function(datW, truth) {
  out <- list()

  blTrueContrs <- NA
  blTrueMeans <- NA
  BVTrueMeans <- NA # avoid NOTE in R CMD check
  e <- list2env(truth, envir = environment())
  levSpecies <- names(blTrueContrs)

  blObsMeans <- list()
  for (species in levSpecies) {
    tmp <- matrix(NA, nrow(blTrueMeans[[species]]), ncol(blTrueMeans[[species]]),
      dimnames = dimnames(blTrueMeans[[species]])
    )
    for (block in rownames(tmp)) {
      for (type in colnames(tmp)) {
        idx <- which(datW$block == block & datW$type == type)
        tmp[block, type] <- mean(datW[idx, paste0("yield_", species)], na.rm = TRUE)
      }
    }
    blObsMeans[[species]] <- tmp
  }
  out[["blObsMeans"]] <- blObsMeans

  BVObsMeans <- list()
  for (species in levSpecies) {
    otherSpecies <- levSpecies[levSpecies != species]
    tmp <- array(
      dim = dim(BVTrueMeans[[species]]),
      dimnames = dimnames(BVTrueMeans[[species]])
    )
    for (geno in dimnames(tmp)[[1]]) {
      for (type in dimnames(tmp)[[2]]) {
        tBV <- "DBV"
        idx <- which(datW[[paste0("geno_", species)]] == geno &
          datW$type == type)
        tmp[geno, type, tBV] <- mean(datW[idx, paste0("yield_", species)], na.rm = TRUE)
        tBV <- "SBV"
        idx <- which(datW[[paste0("geno_", species)]] == geno &
          datW$type == type)
        tmp[geno, type, tBV] <- mean(datW[idx, paste0("yield_", otherSpecies)], na.rm = TRUE)
      }
    }
    BVObsMeans[[species]] <- tmp
  }
  out[["BVObsMeans"]] <- BVObsMeans

  idxIC <- which(!is.na(datW$geno_S1) & !is.na(datW$geno_S2))
  blObsContrs <- list()
  for (species in levSpecies) {
    tmp <- matrix(
      nrow = nlevels(datW$block), ncol = 2,
      dimnames = list(
        c("(Intercept)", paste0("block", levels(datW$block)[-1])),
        c("SC", "IC")
      )
    )
    for (type in c("SC", "IC")) {
      tmp[, type] <- solve(cbind(1, contr.sum(levels(datW[idxIC, "block"])))) %*%
        +blObsMeans[[species]][, type]
    }
    blObsContrs[[species]] <- tmp
  }
  out[["blObsContrs"]] <- blObsContrs

  BVObsContrs <- list()
  for (species in levSpecies) {
    otherSpecies <- levSpecies[levSpecies != species]
    tmp <- array(
      dim = dim(BVObsMeans[[species]]),
      dimnames = list(
        c(
          "(Intercept)",
          paste0(
            "geno_", species,
            levels(datW[[paste0("geno_", species)]])[-1]
          )
        ),
        dimnames(BVObsMeans[[species]])[[2]],
        dimnames(BVObsMeans[[species]])[[3]]
      )
    )
    for (type in c("SC", "IC")) {
      for (tBV in c("DBV", "SBV")) {
        tmp[, type, tBV] <-
          solve(cbind(1, contr.sum(levels(datW[idxIC, paste0("geno_", species)])))) %*%
          +BVObsMeans[[species]][, type, tBV]
      }
    }
    BVObsContrs[[species]] <- tmp
  }
  out[["BVObsContrs"]] <- BVObsContrs

  return(out)
}

##' @noRd
computeYieldICMatricesInterMixDesign <- function(datW, truth, sep) {
  blTrueContrs <- NA
  BV <- NA
  DBVxSBV <- NA
  allErrors <- NA
  B_IC <- NA # avoid NOTE in R CMD check
  e <- list2env(truth, envir = environment())
  levSpecies <- names(blTrueContrs)

  idxIC <- which(!is.na(datW$geno_S1) & !is.na(datW$geno_S2))
  datW_IC <- droplevels(datW[idxIC, ])

  form <- "~ 1"
  lContr <- list()
  if (nlevels(datW_IC$block) > 1) {
    form <- paste0(form, " + block")
    lContr[["block"]] <- "contr.sum"
  }
  form <- paste0(form, " + geno_S2")
  lContr[["geno_S2"]] <- "contr.sum"
  X_IC <- model.matrix(as.formula(form), datW_IC,
    contrasts.arg = lContr
  )

  BV_IC <- cbind(
    "S1" = as.vector(BV$S1[levels(datW_IC$geno_S1), "DBV"]),
    "S2" = as.vector(BV$S1[levels(datW_IC$geno_S1), "SBV"])
  )
  rownames(BV_IC) <- levels(datW_IC$geno_S1)
  Z_DS_IC <- model.matrix(~ 0 + geno_S1, datW_IC)
  colnames(Z_DS_IC) <- gsub("^geno_S1", "", colnames(Z_DS_IC))
  stopifnot(all(colnames(Z_DS_IC) == rownames(BV_IC)))

  datW_IC$mixID <- factor(paste0(datW_IC$geno_S2, sep, datW_IC$geno_S1))
  Z_DxS_IC <- model.matrix(~ 0 + mixID, datW_IC)
  colnames(Z_DxS_IC) <- gsub("^mixID", "", colnames(Z_DxS_IC))
  DBVxSBV_IC <- DBVxSBV[colnames(Z_DxS_IC), ]

  E_IC <- do.call(rbind, allErrors[idxIC])

  Y_IC <- X_IC %*% B_IC + Z_DS_IC %*% BV_IC + E_IC
  Y_IC <- Y_IC + Z_DxS_IC %*% DBVxSBV_IC

  try(stopifnot(all.equal(Y_IC[, 1],
    datW_IC$yield_S1,
    check.attributes = FALSE
  )))
  try(stopifnot(all.equal(Y_IC[, 2],
    datW_IC$yield_S2,
    check.attributes = FALSE
  )))

  return(Y_IC)
}

##' @noRd
computeYieldSCMatricesInterMixDesign <- function(datW, truth) {
  out <- list()

  BV <- NA
  SIGVs <- NA
  beta_SC_f <- NA
  beta_SC_t <- NA
  allErrors <- NA # avoid NOTE in R CMD check
  e <- list2env(truth, envir = environment())
  levSpecies <- names(BV)

  datWs_SC <- list()
  for (species in levSpecies) {
    otherSpecies <- levSpecies[levSpecies != species]
    idxSC <- which(!is.na(datW[[paste0("geno_", species)]]) &
      is.na(datW[[paste0("geno_", otherSpecies)]]))
    datW_SC <- droplevels(datW[idxSC, c("standID", "type", "block")])

    X_SC <- matrix(nrow = 0, ncol = 0)
    Z_D <- matrix(nrow = 0, ncol = 0)
    Z_SIGV <- matrix(nrow = 0, ncol = 0)

    form <- "~ 1"
    lContr <- NULL
    if (nlevels(datW_SC$block) > 1) {
      form <- paste0(form, " + block")
      lContr[["block"]] <- "contr.sum"
    }
    if (species == "S1") {
      X_SC <- model.matrix(as.formula(form), datW_SC,
        contrasts.arg = lContr
      )
      beta <- beta_SC_f
      Z_D <- model.matrix(~ 0 + standID, datW_SC)
      colnames(Z_D) <- gsub("^standID", "", colnames(Z_D))
      DBV <- BV$S1[levels(datW_SC$standID), "DBV"]
      stopifnot(all(colnames(Z_D) == names(DBV)))
      SIGV <- SIGVs$S1[levels(datW_SC$standID)]
      stopifnot(all(colnames(Z_D) == names(SIGV)))
    } else if (species == "S2") {
      form <- paste0(form, " + standID")
      lContr[["standID"]] <- "contr.sum"
      X_SC <- model.matrix(as.formula(form), datW_SC,
        contrasts.arg = lContr
      )
      beta <- beta_SC_t
      ## TODO: if SIGVs$S2 is included, can it be distinguished from DBV?
      ## => standID included twice to make X_SC
    }

    E_SC <- do.call(c, allErrors[idxSC])

    y_SC <- X_SC %*% beta
    if (species == "S1") {
      y_SC <- y_SC + Z_D %*% DBV + Z_D %*% SIGV
    }
    y_SC <- y_SC + E_SC
    out[[species]] <- y_SC

    datW_SC$true_yield <- datW[idxSC, paste0("true_yield_", species)]
    datW_SC$yield <- datW[idxSC, paste0("yield_", species)]
    try(stopifnot(all.equal(as.vector(y_SC[, 1]),
      datW_SC$yield,
      check.attributes = FALSE
    )))
  }

  return(out)
}

##' @noRd
addTotalYieldInterMixDesign <- function(datW) {
  datW$tot_yield <- NA
  idx <- which(datW$type2 == "sole_S1")
  datW$tot_yield[idx] <- datW$yield_S1[idx]
  idx <- which(datW$type2 == "sole_S2")
  datW$tot_yield[idx] <- datW$yield_S2[idx]
  idx <- which(datW$type2 == "mix")
  datW$tot_yield[idx] <- datW$yield_S1[idx] + datW$yield_S2[idx]
  return(datW)
}

##' @noRd
convDatW2LInterMixDesign <- function(datW, sep) {
  datL <- pivotMixData2Long(datW,
    colC = "standID", sep = sep,
    genos = list(
      "S1" = as.character(datW$geno_S1),
      "S2" = as.character(datW$geno_S2)
    )
  )
  datL$true_focal_yield <- datL$true_yield_S1
  datL$focal_yield <- datL$yield_S1
  idx <- which(grepl("^gS2_", datL$focal))
  datL$true_focal_yield[idx] <- datL$true_yield_S2[idx]
  datL$focal_yield[idx] <- datL$yield_S2[idx]
  datL$focal_species <- "S1"
  idx <- grep("^gS2_", datL$focal)
  datL$focal_species[idx] <- "S2"
  datL$focal_species <- factor(datL$focal_species)
  return(datL)
}

##' @noRd
addPropsInterMixDesign <- function(datL, sowingDensities) {
  stopifnot(
    all(c("type2", "focal") %in% colnames(datL)),
    is.list(sowingDensities),
    length(sowingDensities) == 2,
    all(names(sowingDensities) == c("S1", "S2")),
    all(sapply(sowingDensities, function(x) {
      all(names(x) == c("SC", "IC"))
    }))
  )

  out <- list()

  overallDensity <- sum(sapply(sowingDensities, `[`, "IC"))
  props <- sapply(sowingDensities, function(x) {
    as.numeric(x["IC"]) / overallDensity
  })
  stopifnot(all.equal(sum(props), 1))
  out$props <- props

  datL$prop <- 1
  idx <- which(datL$type2 == "mix" & grepl("^gS1_", datL$focal))
  stopifnot(length(idx) > 0)
  datL$prop[idx] <- props["S1"]
  idx <- which(datL$type2 == "mix" & grepl("^gS2_", datL$focal))
  stopifnot(length(idx) > 0)
  datL$prop[idx] <- props["S2"]
  out$datL <- datL

  return(out)
}

##' Simulate w.r.t. a DBV-SBV model for interspecific mixtures
##'
##' Simulates an incomplete (sparse) yet balanced, tester-based design with respect to a DBV-SBV model for interspecific mixtures.
##' @param GRMs named list with a genomic relationship matrix per species, named S1 and S2; should be a diagonal for S2 if it is a tester (see \code{parsDesign})
##' @param parsGenetics list of genetic parameters
##' @param parsDesign list of design parameters
##' @param sep character separating the names of a mixture components
##' @param seed seed for the generation of pseudo-random numbers
##' @return list
##' @seealso \code{\link{fitDBVSBVinter}}
##' @author Timothee Flutre
##' @examples
##' ## simulate a data set with both sole crops and intercrops:
##' GRMs <- list("S1" = diag(100),
##'              "S2" = diag(2))
##' dimnames(GRMs$S1) <- list(paste0("gS1_", 1:100), paste0("gS1_", 1:100))
##' dimnames(GRMs$S2) <- list(paste0("gS2_", 1:2), paste0("gS2_", 1:2))
##' out <- simulDBVSBVinter(GRMs)
##' names(out)
##' str(out$datW)
##' str(out$datL)
##'
##' ## see the third vignette for more details
##' @export
simulDBVSBVinter <- function(GRMs,
                             parsGenetics = list(
                               mu = list(
                                 "S1" = c("SC" = 65, "IC" = 32),
                                 "S2" = c("SC" = 30, "IC" = 27)
                               ),
                               sdBlocks = 4,
                               CV_g = c("S1" = 0.08, "S2" = 0.08),
                               prop_var_SBV = c("S1" = 0.2, "S2" = 0.2),
                               prop_var_SIGV = c("S1" = 0.5, "S2" = 0),
                               cor_DBV_SBV = c("S1" = -0.9, "S2" = -0.9),
                               prop_var_DBVxSBV = c("S1" = 0, "S2" = 0),
                               use_GRM_for_DBVxSBV = TRUE,
                               cor_DxS = 0,
                               H2_SC = c("S1" = 0.7, "S2" = 0.7),
                               T2 = c("S1" = 0.7, "S2" = 0.7),
                               cor_err_IC = -0.2
                             ),
                             parsDesign = list(
                               nbBlocks = 2,
                               nbPlotsPerBl = 500,
                               nbYs = 10,
                               sowingDensities = list(
                                 "S1" = c("SC" = 300, "IC" = 150),
                                 "S2" = c("SC" = 40, "IC" = 40)
                               ),
                               testersS2 = TRUE,
                               incomplete_balanced = TRUE
                             ),
                             sep = "+",
                             seed = NULL) {
  e <- list2env(parsDesign, envir = environment())

  stopifnot(
    is.list(GRMs),
    length(GRMs) == 2,
    all(sapply(GRMs, function(GRM) {
      nrow(GRM) == ncol(GRM)
    })),
    !is.null(names(GRMs)),
    all(c("S1", "S2") == names(GRMs)),
    all(sapply(GRMs, function(GRM) {
      !is.null(dimnames(GRM))
    })),
    parsDesign$nbBlocks <= 26
  )
  if (parsDesign$testersS2) {
    stopifnot(all(diag(GRMs$S2) == 1))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  nbGenos <- sapply(GRMs, nrow)
  levGenos <- lapply(GRMs, rownames)
  levBlocks <- LETTERS[1:parsDesign$nbBlocks]

  dat_bl1 <- makeCompleteDesignOneBlockInterMixDesign(levGenos, sep)
  datW <- makeCompleteDesignInterMixDesign(dat_bl1, levBlocks)
  if (parsDesign$incomplete_balanced) {
    datW <- allocGenosMixsInterMixDesign(datW, parsDesign$nbBlocks, levGenos)
  }
  datW <- setMplotCoordsInterMixDesign(datW, parsDesign$nbYs)
  truth <- drawParamsInterMixDesign(datW, GRMs, sep, parsGenetics, parsDesign)
  datW <- computeYieldForLoopInterMixDesign(datW, truth, sep)
  Y_IC <- computeYieldICMatricesInterMixDesign(datW, truth, sep)
  y_SC <- computeYieldSCMatricesInterMixDesign(datW, truth)
  obsMC <- computeObsMeansContrsInterMixDesign(datW, truth)
  datW <- addTotalYieldInterMixDesign(datW)
  datL <- convDatW2LInterMixDesign(datW, sep)
  tmp <- addPropsInterMixDesign(datL, parsDesign$sowingDensities)
  datL <- tmp$datL
  props <- tmp$props

  return(list(
    truth = truth, datW = datW, datL = datL, obsMC = obsMC,
    sowingDensities = parsDesign$sowingDensities, props = props
  ))
}
