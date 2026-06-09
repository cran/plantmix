fitLM_prepIn <- function(y, X) {
  out <- list()

  listData <- list(
    model = "lm",
    y = y,
    X = X
  )
  out$listData <- listData

  listParams <- list(
    beta = rep(0, ncol(X)),
    log_sigma = log(1)
  )
  out$listParams <- listParams

  out$vecRnd <- NULL

  return(out)
}

fitLM_prepOut <- function(outTmb) {
  out <- list()

  sdrep <- outTmb$sdrep
  p <- length(sdrep$par.fixed) - 1

  out$beta <- data.frame(
    estim = sdrep$par.fixed[1:p],
    stderr = sqrt(diag(sdrep$cov.fixed)[1:p])
  )
  rownames(out$beta) <- NULL

  ## Delta method: Var[X] ~= Var[log(X)] E[X]^2
  ## https://stats.stackexchange.com/a/488166/3459
  out$sigma <- data.frame(
    estim = exp(sdrep$par.fixed[p + 1]),
    stderr = diag(sdrep$cov.fixed)[p + 1]
  )
  out$sigma$stderr <- sqrt(out$sigma$stderr * out$sigma$estim^2)
  rownames(out$sigma) <- NULL

  return(out)
}

##' Linear model
##'
##' Fits a linear model.
##' @param y vector of responses
##' @param X design matrix
##' @param verbose verbosity level
##' @return list of parameter estimates
##' @author Timothee Flutre
##' @examples
##' \dontrun{## fake data
##' set.seed(12345)
##' n <- 100
##' p <- 5
##' X <- matrix(rnorm(n*p), n, p)
##' (beta <- rnorm(p, sd=2))
##' sigma <- 1
##' y <- rnorm(n=n, mean=X %*% beta, sd=sigma)
##'
##' ## model fit:
##' (fit <- fitLM(y, X))
##'
##' ## comparison:
##' fit_ref <- lm(y ~ 0 + X, data.frame(y=y, X))
##' coef(summary(fit_ref))
##' sigma(fit_ref)
##' }
##' @noRd
fitLM <- function(y, X, verbose = 0) {
  ## quick reformat and check
  if (is(y, "matrix")) {
    if (ncol(y) == 1) {
      y <- as.vector(y)
    }
  }
  stopifnot(
    is.vector(y),
    is.matrix(X)
  )

  out <- list()

  ## input preparation
  inputs4TMB <- fitLM_prepIn(y, X)

  ## automatic differentiation
  f <- MakeADFun(
    data = inputs4TMB$listData,
    parameters = inputs4TMB$listParams,
    random = inputs4TMB$vecRnd,
    DLL = "plantmix_TMBExports",
    silent = ifelse(verbose == 0, TRUE, FALSE)
  )
  out$f <- f

  ## optimization
  fit <- nlminb(start = f$par, objective = f$fn, gradient = f$gr, hessian = NULL)
  out$fit <- fit

  ## uncertainty quantification
  sdrep <- sdreport(f)
  out$sdrep <- sdrep

  ## output preparation
  out <- fitLM_prepOut(out)

  return(out)
}
