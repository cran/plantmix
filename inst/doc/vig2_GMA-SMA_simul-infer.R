## ----echo=FALSE---------------------------------------------------------------
suppressPackageStartupMessages(library(knitr))
opts_chunk$set(
  echo = TRUE, warning = TRUE, message = TRUE, cache = FALSE,
  fig.align = "center", collapse = TRUE
)
opts_knit$set(progress = TRUE, verbose = TRUE)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(plantmix))
suppressPackageStartupMessages(library(lme4))

## ----echo=FALSE---------------------------------------------------------------
## parallelize to speed-up usage of MM4LMM
## suppressPackageStartupMessages(library(parallel))
## nbCores <- detectCores() - 1

## ----time_0, echo=FALSE-------------------------------------------------------
## Execution time (see the appendix):
t0 <- proc.time()

## -----------------------------------------------------------------------------
set.seed(12345)

## -----------------------------------------------------------------------------
nbGenos <- 200
levGenos <- sprintf(
  fmt = paste0("g%0", floor(log10(nbGenos)) + 1, "i"),
  1:nbGenos
)

nbSnps <- 1000
levSnps <- sprintf(
  fmt = paste0("s%0", floor(log10(nbSnps)) + 1, "i"),
  1:nbSnps
)

## -----------------------------------------------------------------------------
nb_pops <- 10
weak_div_pops <- diag(nb_pops)
weak_div_pops[upper.tri(weak_div_pops)] <- 0.9
weak_div_pops[lower.tri(weak_div_pops)] <- weak_div_pops[upper.tri(weak_div_pops)]
tmp <- rep(nbGenos / nb_pops, nb_pops - 1)
tmp <- c(tmp, nbGenos - sum(tmp))
snpGenos <- simulGenosDoseStruct(
  nb_genos = tmp,
  nb_snps = nbSnps,
  div_pops = weak_div_pops,
  geno_IDs = levGenos,
  snp_IDs = levSnps
)
dim(snpGenos)
snpGenos[1:3, 1:4]

## -----------------------------------------------------------------------------
GRM <- estimGRM(snpGenos)
GRM <- as.matrix(Matrix::nearPD(GRM)$mat)

## ----class.source='fold-hide'-------------------------------------------------
image(Matrix(GRM), main = "GRM")
hist(diag(GRM), main = "GRM")
hist(GRM[upper.tri(GRM)], main = "GRM")

## -----------------------------------------------------------------------------
nbMixes <- 1000
design <- getDesignBinaryVarMix(levGenos, nbMixes = nbMixes)
tmp <- getMixturesPerGeno(getMixtureList(design$combs))
table(sapply(tmp, length)) # each genotype is in the same nb of mixtures -> balanced design

## ----fig.width=8, echo=FALSE--------------------------------------------------
## plotDesignVarMix(design$graph, levGenos, subplots = "diallel")
plotDiallel(design$diallel, main = paste0(
  nbGenos, " genotypes and ",
  nbMixes, " binary mixtures"
))

## -----------------------------------------------------------------------------
monoStands <- paste(levGenos, levGenos, sep = "_")

pairs <- design$combs
mixStands <- paste(pairs[, 1], pairs[, 2], sep = "_")
IDs <- c(monoStands, mixStands)

nbBlocks <- 3
blocks <- LETTERS[1:nbBlocks]

dat <- do.call(rbind, lapply(blocks, function(block) {
  data.frame(
    ID = IDs,
    block = block,
    stringsAsFactors = TRUE
  )
}))

## -----------------------------------------------------------------------------
listContr <- list(block = "contr.sum")
X <- model.matrix(~ 1 + block, data = dat, contrasts = listContr)

allZ <- list()
allZ$GMA <- mkZGMA(dat, col = "ID", sep = "_")
system.time(
  allZ <- append(
    allZ,
    mkAllZSMA(dat, col = "ID", sep = "_", verbose = 1)
  )
)
sapply(allZ, dim)

## -----------------------------------------------------------------------------
truth <- list(
  "intercept" = 100,
  "var_GMA" = 10,
  "var_SMA" = 2,
  "var_SMA_ij" = 1.5,
  "var_SMA_ii" = 0.8,
  "var_error" = 1
)
set.seed(1234)
truth[["blockEffs"]] <- sample(x = c(-1, 1), size = nbBlocks - 1, replace = TRUE) *
  rnorm(n = nbBlocks - 1, mean = 3, sd = 2)
GRM <- diag(nbGenos)
truth[["GMAs"]] <- MASS::mvrnorm(
  n = 1, mu = rep(0, nbGenos),
  Sigma = truth$var_GMA * GRM
)
truth[["SMAs"]] <- rnorm(n = length(IDs), mean = 0, sd = sqrt(truth$var_SMA))
truth[["SMAs_ij"]] <- rnorm(n = length(mixStands), mean = 0, sd = sqrt(truth$var_SMA_ij))
truth[["SMAs_ii"]] <- rnorm(n = length(monoStands), mean = 0, sd = sqrt(truth$var_SMA_ii))
truth[["errors"]] <- rnorm(n = nrow(dat), mean = 0, sd = sqrt(truth$var_error))

## -----------------------------------------------------------------------------
y1 <- X %*% c(truth$intercept, truth$blockEffs) +
  allZ$GMA %*% truth$GMAs +
  truth$errors
dat$pheno1 <- y1[, 1]

y2 <- X %*% c(truth$intercept, truth$blockEffs) +
  allZ$GMA %*% truth$GMAs +
  allZ$SMA_mod2 %*% truth$SMAs +
  truth$errors
dat$pheno2 <- y2[, 1]

y3 <- X %*% c(truth$intercept, truth$blockEffs) +
  allZ$GMA %*% truth$GMAs +
  allZ$SMA_mod3 %*% truth$SMAs +
  truth$errors
dat$pheno3 <- y3[, 1]

y2p <- X %*% c(truth$intercept, truth$blockEffs) +
  allZ$GMA %*% truth$GMAs +
  allZ$SMA_mod2p %*% truth$SMAs_ij +
  truth$errors
dat$pheno2p <- y2p[, 1]

y2pp <- X %*% c(truth$intercept, truth$blockEffs) +
  allZ$GMA %*% truth$GMAs +
  allZ$SMA_mod2pp_ij %*% truth$SMAs_ij +
  allZ$SMA_mod2pp_ii %*% truth$SMAs_ii +
  truth$errors
dat$pheno2pp <- y2pp[, 1]

y3p <- X %*% c(truth$intercept, truth$blockEffs) +
  allZ$GMA %*% truth$GMAs +
  allZ$SMA_mod3p_ij %*% truth$SMAs_ij +
  allZ$SMA_mod3p_ii %*% truth$SMAs_ii +
  truth$errors
dat$pheno3p <- y3p[, 1]

## -----------------------------------------------------------------------------
pkgs <- c("lme4", "TMB") # , "MM4LMM")
if ("MM4LMM" %in% pkgs) {
  suppressPackageStartupMessages(library(MM4LMM))
}

runAllPkgs <- function(pkgs, form, dat, listZ, listVCov, listContr, ...) {
  ## mcMap(function(pkg) {
  Map(function(pkg) {
    print(paste0("fit with ", pkg, "..."))
    fitGMASMA(form, dat, listZ, pkg, listVCov, listContr)
    st <- system.time(
      try(fit <- fitGMASMA(form, dat, listZ, pkg, listVCov, listContr, REML = TRUE, ...))
    )
    print(st)
    fit
  }, pkgs) # , mc.cores = nbCores)
}

## -----------------------------------------------------------------------------
form <- pheno1 ~ 1 + block
listZ <- list("GMA" = allZ$GMA)
listVCov <- list("GMA" = GRM)
listContr <- list(block = "contr.sum")
fits1 <- runAllPkgs(pkgs, form, dat, listZ, listVCov, listContr)

## -----------------------------------------------------------------------------
print("lme4:")
tmp <- as.data.frame(VarCorr(fits1$lme4))
tmp <- setNames(tmp$vcov, tmp$grp)
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_error")]),
  estim = tmp[c("GMA", "Residual")]
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

print("TMB:")
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_error")]),
  estim = do.call(c, fits1$TMB$report[c("var_GMA", "var_err")])
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

if ("MM4LMM" %in% names(fits1)) {
  print(fits1$MM4LMM$Sigma2)
}

## ----echo=FALSE---------------------------------------------------------------
if ("INLA" %in% names(fits1)) {
  m <- fits1$INLA$internal.marginals.hyperpar[[1]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_error
  m <- fits1$INLA$internal.marginals.hyperpar[[2]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_GMA
}

## -----------------------------------------------------------------------------
form <- pheno2 ~ 1 + block
listZ <- list(
  "GMA" = allZ$GMA,
  "SMA" = allZ$SMA_mod2
)
listVCov <- list(
  "GMA" = GRM,
  "SMA" = diag(ncol(listZ$SMA))
)
listContr <- list(block = "contr.sum")
fits2 <- runAllPkgs(pkgs, form, dat, listZ, listVCov, listContr)

## -----------------------------------------------------------------------------
print("lme4:")
tmp <- as.data.frame(VarCorr(fits2$lme4))
tmp <- setNames(tmp$vcov, tmp$grp)
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA", "var_error")]),
  estim = tmp[c("GMA", "SMA", "Residual")]
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

print("TMB:")
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA", "var_error")]),
  estim = do.call(c, fits2$TMB$report[c("var_GMA", "var_SMA1", "var_err")])
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

if ("MM4LMM" %in% names(fits2)) {
  print(fits2$MM4LMM$Sigma2)
}

## ----echo=FALSE---------------------------------------------------------------
if ("INLA" %in% names(fits2)) {
  ## fits2$INLA$summary.hyperpar
  ## names(fits2$INLA$internal.marginals.hyperpar)
  m <- fits2$INLA$internal.marginals.hyperpar[[1]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_error
  m <- fits2$INLA$internal.marginals.hyperpar[[2]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_GMA
  m <- fits2$INLA$internal.marginals.hyperpar[[3]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_GMA
}

## -----------------------------------------------------------------------------
form <- pheno3 ~ 1 + block
listZ <- list(
  "GMA" = allZ$GMA,
  "SMA" = allZ$SMA_mod3
)
listVCov <- list(
  "GMA" = GRM,
  "SMA" = diag(ncol(listZ$SMA))
)
listContr <- list(block = "contr.sum")
fits3 <- runAllPkgs(pkgs, form, dat, listZ, listVCov, listContr)

## -----------------------------------------------------------------------------
print("lme4:")
tmp <- as.data.frame(VarCorr(fits3$lme4))
tmp <- setNames(tmp$vcov, tmp$grp)
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA", "var_error")]),
  estim = tmp[c("GMA", "SMA", "Residual")]
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

print("TMB:")
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA", "var_error")]),
  estim = do.call(c, fits3$TMB$report[c("var_GMA", "var_SMA1", "var_err")])
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

if ("MM4LMM" %in% names(fits3)) {
  print(fits3$MM4LMM$Sigma2)
}

## ----echo=FALSE---------------------------------------------------------------
if ("INLA" %in% names(fits3)) {
  ## fits3$INLA$summary.hyperpar
  ## names(fits3$INLA$internal.marginals.hyperpar)
  m <- fits3$INLA$internal.marginals.hyperpar[[1]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_error
  m <- fits3$INLA$internal.marginals.hyperpar[[2]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_GMA
}

## -----------------------------------------------------------------------------
form <- pheno2p ~ 1 + block
listZ <- list(
  "GMA" = allZ$GMA,
  "SMA" = allZ$SMA_mod2p
)
listVCov <- list(
  "GMA" = GRM,
  "SMA" = diag(ncol(listZ$SMA))
)
listContr <- list(block = "contr.sum")
fits2p <- runAllPkgs(pkgs, form, dat, listZ, listVCov, listContr)

## -----------------------------------------------------------------------------
print("lme4:")
tmp <- as.data.frame(VarCorr(fits2p$lme4))
tmp <- setNames(tmp$vcov, tmp$grp)
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA", "var_error")]),
  estim = tmp[c("GMA", "SMA", "Residual")]
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

print("TMB:")
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA", "var_error")]),
  estim = do.call(c, fits2p$TMB$report[c("var_GMA", "var_SMA1", "var_err")])
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

if ("MM4LMM" %in% names(fits2p)) {
  print(fits2p$MM4LMM$Sigma2)
}

## ----echo=FALSE---------------------------------------------------------------
if ("INLA" %in% names(fits2p)) {
  ## fits2p$INLA$summary.hyperpar
  ## names(fits2p$INLA$internal.marginals.hyperpar)
  m <- fits2p$INLA$internal.marginals.hyperpar[[1]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_error
  m <- fits2p$INLA$internal.marginals.hyperpar[[2]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_GMA
  m <- fits2p$INLA$internal.marginals.hyperpar[[3]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_GMA
}

## -----------------------------------------------------------------------------
form <- pheno2pp ~ 1 + block
listZ <- list(
  "GMA" = allZ$GMA,
  "SMA_ij" = allZ$SMA_mod2pp_ij,
  "SMA_ii" = allZ$SMA_mod2pp_ii
)
listVCov <- list(
  "GMA" = GRM,
  "SMA_ij" = diag(ncol(listZ$SMA_ij)),
  "SMA_ii" = diag(ncol(listZ$SMA_ii))
)
listContr <- list(block = "contr.sum")
fits2pp <- runAllPkgs(pkgs, form, dat, listZ, listVCov, listContr)

## -----------------------------------------------------------------------------
print("lme4:")
tmp <- as.data.frame(VarCorr(fits2pp$lme4))
tmp <- setNames(tmp$vcov, tmp$grp)
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA_ij", "var_SMA_ii", "var_error")]),
  estim = tmp[c("GMA", "SMA_ij", "SMA_ii", "Residual")]
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

print("TMB:")
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA_ij", "var_SMA_ii", "var_error")]),
  estim = do.call(c, fits2pp$TMB$report[c("var_GMA", "var_SMA1", "var_SMA2", "var_err")])
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

if ("MM4LMM" %in% names(fits2pp)) {
  print(fits2pp$MM4LMM$Sigma2)
}

## ----echo=FALSE---------------------------------------------------------------
if ("INLA" %in% names(fits2pp)) {
  ## fits2pp$INLA$summary.hyperpar
  ## names(fits2pp$INLA$internal.marginals.hyperpar)
  m <- fits2pp$INLA$internal.marginals.hyperpar[[1]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_error
  m <- fits2pp$INLA$internal.marginals.hyperpar[[2]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_GMA
  m <- fits2pp$INLA$internal.marginals.hyperpar[[3]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_SMA_ij
  m <- fits2pp$INLA$internal.marginals.hyperpar[[4]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_SMA_ii
}

## -----------------------------------------------------------------------------
form <- pheno3p ~ 1 + block
listZ <- list(
  "GMA" = allZ$GMA,
  "SMA_ij" = allZ$SMA_mod3p_ij,
  "SMA_ii" = allZ$SMA_mod3p_ii
)
listVCov <- list(
  "GMA" = GRM,
  "SMA_ij" = diag(ncol(listZ$SMA_ij)),
  "SMA_ii" = diag(ncol(listZ$SMA_ii))
)
listContr <- list(block = "contr.sum")
fits3p <- runAllPkgs(pkgs, form, dat, listZ, listVCov, listContr)

## -----------------------------------------------------------------------------
print("lme4:")
tmp <- as.data.frame(VarCorr(fits3p$lme4))
tmp <- setNames(tmp$vcov, tmp$grp)
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA_ij", "var_SMA_ii", "var_error")]),
  estim = tmp[c("GMA", "SMA_ij", "SMA_ii", "Residual")]
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

print("TMB:")
checks <- data.frame(
  truth = do.call(c, truth[c("var_GMA", "var_SMA_ij", "var_SMA_ii", "var_error")]),
  estim = do.call(c, fits3p$TMB$report[c("var_GMA", "var_SMA1", "var_SMA2", "var_err")])
)
checks$nBE <- normBiasError(checks$estim, checks$truth)
print(checks)

if ("MM4LMM" %in% names(fits3p)) {
  print(fits3p$MM4LMM$Sigma2)
}

## ----echo=FALSE---------------------------------------------------------------
if ("INLA" %in% names(fits3p)) {
  ## fits3p$INLA$summary.hyperpar
  ## names(fits3p$INLA$internal.marginals.hyperpar)
  m <- fits3p$INLA$internal.marginals.hyperpar[[1]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_error
  m <- fits3p$INLA$internal.marginals.hyperpar[[2]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_GMA
  m <- fits3p$INLA$internal.marginals.hyperpar[[3]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_SMA_ij
  m <- fits3p$INLA$internal.marginals.hyperpar[[4]]
  m.var <- INLA::inla.tmarginal(function(x) 1 / exp(x), m)
  INLA::inla.zmarginal(m.var) # var_SMA_ii
}

## -----------------------------------------------------------------------------
t1 <- proc.time()
t1 - t0
print(sessionInfo(), locale = FALSE)

