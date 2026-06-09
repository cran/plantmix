## ----echo=FALSE---------------------------------------------------------------
suppressPackageStartupMessages(library(knitr))
opts_chunk$set(
  echo = TRUE, warning = TRUE, message = TRUE, cache = FALSE,
  fig.align = "center", collapse = TRUE
)
opts_knit$set(progress = TRUE, verbose = TRUE)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(plantmix))

## ----time_0, echo=FALSE-------------------------------------------------------
## Execution time (see the appendix):
t0 <- proc.time()

## ----echo=FALSE---------------------------------------------------------------
dat <- data.frame(
  id = c("stand1", "stand2", "stand3"),
  genos = c("g1", "g1-g2", "g1-g2-g3"),
  type = c("monovarietal", "mixed", "mixed"),
  order = c(1, 2, 3)
)
knitr::kable(dat, row.names = TRUE)

## ----echo=FALSE---------------------------------------------------------------
Z_G <- mkZGMA(df = dat, col = "genos", sep = "-")
rownames(Z_G) <- dat$genos
knitr::kable(as.data.frame(Z_G), digits = 2, row.names = TRUE)

## ----echo=FALSE---------------------------------------------------------------
Z_S <- mkZSMA(df = dat, col = "genos", sep = "-", inc_SMA_ii = "only_pur")
rownames(Z_S) <- dat$genos
knitr::kable(as.data.frame(Z_S), digits = 2, row.names = TRUE)

## ----echo=FALSE---------------------------------------------------------------
Z_S <- mkZSMA(df = dat, col = "genos", sep = "-", inc_SMA_ii = "pur_mix")
rownames(Z_S) <- dat$genos
knitr::kable(as.data.frame(Z_S), digits = 2, row.names = TRUE)

## ----echo=FALSE---------------------------------------------------------------
Z_S <- mkZSMA(df = dat, col = "genos", sep = "-", inc_SMA_ii = "no")
rownames(Z_S) <- dat$genos
knitr::kable(as.data.frame(Z_S), digits = 2, row.names = TRUE)

## ----echo=FALSE---------------------------------------------------------------
Z_S <- mkZSMA(df = dat, col = "genos", sep = "-", inc_SMA_ii = "only_pur")
rownames(Z_S) <- dat$genos
isMono <- (sapply(strsplit(colnames(Z_S), "-"), anyDuplicated) == 2)

## ----echo=FALSE---------------------------------------------------------------
Z_S_ij <- Z_S[, !isMono]
knitr::kable(as.data.frame(Z_S_ij), digits = 2, row.names = TRUE)

## ----echo=FALSE---------------------------------------------------------------
Z_S_ii <- Z_S[, isMono]
knitr::kable(as.data.frame(Z_S_ii), digits = 2, row.names = TRUE)

## ----echo=FALSE---------------------------------------------------------------
Z_S <- mkZSMA(df = dat, col = "genos", sep = "-", inc_SMA_ii = "pur_mix")
rownames(Z_S) <- dat$genos
isMono <- (sapply(strsplit(colnames(Z_S), "-"), anyDuplicated) == 2)

## ----echo=FALSE---------------------------------------------------------------
Z_S_ij <- Z_S[, !isMono]
knitr::kable(as.data.frame(Z_S_ij), digits = 2, row.names = TRUE)

## ----echo=FALSE---------------------------------------------------------------
Z_S_ii <- Z_S[, isMono]
knitr::kable(as.data.frame(Z_S_ii), digits = 2, row.names = TRUE)

## -----------------------------------------------------------------------------
nbGenos <- 25
levGenos <- sprintf(
  fmt = paste0("geno%0", floor(log10(nbGenos)) + 1, "i"),
  1:nbGenos
)
nbMixes <- 75 # only binary and balanced
design <- getDesignBinaryVarMix(levGenos, nbMixes, seed = 12345)
tmp <- getMixturesPerGeno(getMixtureList(design$combs))
table(sapply(tmp, length)) # each genotype is in the same nb of mixtures -> balanced design

## ----echo=FALSE---------------------------------------------------------------
## plotDesignVarMix(design$graph, levGenos)
plotDiallel(design$diallel, main = paste0(
  nbGenos, " genotypes and ",
  nbMixes, " binary mixtures"
))

## -----------------------------------------------------------------------------
nbBlocks <- 3
levBlocks <- LETTERS[1:nbBlocks]
dat <- do.call(rbind, lapply(levBlocks, function(block) {
  data.frame(
    stand = paste0(design$combs$comp1, "_", design$combs$comp2),
    block = block,
    stringsAsFactors = TRUE
  )
}))
str(dat)

## -----------------------------------------------------------------------------
listContr <- list(block = "contr.sum")
X <- model.matrix(~ 1 + block, data = dat, contrasts = listContr)
Z_GMA <- mkZGMA(dat, "stand", sep = "_")
Z_SMA <- mkZSMA(dat, "stand", sep = "_", inc_SMA_ii = "no")

## -----------------------------------------------------------------------------
truth <- list(
  "intercept" = 80,
  "var_GMA" = 10,
  "var_SMA" = 2,
  "var_error" = 1
)
set.seed(1234)
truth[["blockEffs"]] <- sample(x = c(-1, 1), size = nbBlocks - 1, replace = TRUE) *
  rnorm(n = nbBlocks - 1, mean = 3, sd = 5)
truth[["GMA"]] <- rnorm(n = nbGenos, mean = 0, sd = sqrt(truth$var_GMA))
truth[["SMA"]] <- rnorm(n = nbMixes, mean = 0, sd = sqrt(truth$var_SMA))
truth[["errors"]] <- rnorm(n = nrow(dat), mean = 0, sd = sqrt(truth$var_error))

## -----------------------------------------------------------------------------
y <- X %*% c(truth$intercept, truth$blockEffs) +
  Z_GMA %*% truth$GMA + Z_SMA %*% truth$SMA +
  truth$errors
dat$yield <- y[, 1]

boxplot(yield ~ block, data = dat, las = 1, main = "Simulated data")

## -----------------------------------------------------------------------------
listZ1 <- list("GMA" = mkZGMA(dat, "stand", sep = "_"))
system.time(
  fit1 <- fitGMASMA(yield ~ 1 + block, dat, listZ1, pkg = "lme4", contrasts = listContr)
)

listZ2 <- list(
  "GMA" = mkZGMA(dat, "stand", sep = "_"),
  "SMA" = mkZSMA(dat, "stand", sep = "_", inc_SMA_ii = "no")
)
system.time(
  fit2 <- fitGMASMA(yield ~ 1 + block, dat, listZ2, pkg = "lme4", contrasts = , listContr)
)

## ----fig.width=12-------------------------------------------------------------
BLUEs <- fixef(fit2)
data.frame(
  "true" = c(truth$intercept, truth$blockEffs),
  "estim" = BLUEs
)

estV <- as.data.frame(lme4::VarCorr(fit2))
data.frame(
  "true" = c(truth$var_GMA, truth$var_SMA, truth$var_error),
  "estim" = c(estV$vcov[estV$grp == "GMA"], estV$vcov[estV$grp == "SMA"], estV$vcov[estV$grp == "Residual"]),
  row.names = c("GMA", "SMA", "error")
)

BLUPs <- ranef(fit2)
cor(truth$GMA, BLUPs$GMA[, "(Intercept)"])
cor(truth$SMA, BLUPs$SMA[, "(Intercept)"])
op <- par(mfrow = c(1, 2))
for (MA in c("GMA", "SMA")) {
  plot(BLUPs[[MA]][, "(Intercept)"], truth[[MA]],
    xlab = paste0("BLUP(", MA, ")"), ylab = paste0("true ", MA),
    main = "Accuracy with lme4", las = 1, pch = 19
  )
  abline(a = 0, b = 1, v = 0, h = 0, lty = 2)
  abline(lm(truth[[MA]] ~ BLUPs[[MA]][, "(Intercept)"]), col = "red")
}
par(op)

## -----------------------------------------------------------------------------
fits <- list("mod1" = fit1, "mod2" = fit2)
t(summarizeGMASMAs(fits))

## -----------------------------------------------------------------------------
t1 <- proc.time()
t1 - t0
print(sessionInfo(), locale = FALSE)

