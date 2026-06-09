## Functions allowing to run lme4::lmer() with custom design matrices
## of random effects, useful for the statistical analysis of mixtures

##' Linear mixed model
##'
##' Same as R/lme4::lmer() with the addition of \code{myblist}
##' @param formula see \code{lmer}
##' @param data see \code{lmer}; \code{data} and \code{myblist} should be coherent, e.g., be cautious when keeping missing data in the response column of \code{data}
##' @param REML see \code{lmer}
##' @param control see \code{lmer}
##' @param start see \code{lmer}
##' @param verbose see \code{lmer}
##' @param subset see \code{lmer}
##' @param weights see \code{lmer}
##' @param na.action see \code{lmer}
##' @param offset see \code{lmer}
##' @param contrasts see \code{lmer}
##' @param devFunOnly see \code{lmer}
##' @param myblist list which names correspond to random effect(s) for which we want to use a custom design matrix; each component is a list with four components named "ff", "sm", "nl" and "cnms"
##' @return see \code{lmer}
##' @seealso \code{\link{mkZGMA}}
##' @author Maxence Remerand, Timothee Flutre
##' @examples
##' \dontrun{
##' ## generate fake data
##' nbGenos <- 25
##' genos <- sprintf("g%02i", 1:nbGenos)
##' pairs <- t(combn(x=genos, m=2))
##' stands <- paste(pairs[,1], pairs[,2], sep="_")
##' nbBlocks <- 3
##' blocks <- LETTERS[1:nbBlocks]
##' dat <- do.call(rbind, lapply(blocks, function(block){
##'   cbind(stands=as.data.frame(stands, stringsAsFactors=TRUE),
##'         block=as.factor(block))
##' }))
##' X <- model.matrix(~ 1 + block, data=dat, contrasts=list(block="contr.sum"))
##' Z_GMA <- mkZGMA(df=dat, col="stands", sep="_")
##' truth <- list("intercept"=100,
##'               "var_GMA"=10,
##'               "var_error"=1)
##' set.seed(1234)
##' truth[["blockEffs"]] <- sample(x=c(-1,1), size=nbBlocks - 1, replace=TRUE) *
##'   rnorm(n=nbBlocks - 1, mean=3, sd=5)
##' truth[["GMAs"]] <- rnorm(n=nbGenos, mean=0, sd=sqrt(truth$var_GMA))
##' truth[["errors"]] <- rnorm(n=nrow(dat), mean=0, sd=sqrt(truth$var_error))
##' y <- X %*% c(truth$intercept, truth$blockEffs) +
##'   Z_GMA %*% truth$GMAs +
##'   truth$errors
##' dat$pheno <- y[,1]
##' if(FALSE){
##'   hist(dat$pheno, las=1)
##'   boxplot(pheno ~ block, data=dat, las=1)
##' }
##'
##' ## fit the model
##' dat2 <- cbind(dat, "GMA"=dat$stands)
##' myformula <- pheno ~ 1 + block + (1|GMA)
##' myblist <- list(list("ff"=factor(colnames(Z_GMA)),
##'                      "sm"=Matrix::Matrix(t(Z_GMA), sparse=TRUE),
##'                      "nl"=as.integer(ncol(Z_GMA)),
##'                      "cnms"="(Intercept)"))
##' names(myblist) <- lme4:::barnames(lme4::findbars(lme4:::RHSForm(myformula)))
##' fit <- lmerZ(formula=myformula, data=dat2,
##'              REML=TRUE, myblist=myblist)
##'
##' ## check the results
##' as.data.frame(lme4::VarCorr(fit))
##' if(FALSE){
##'   plot(x=lme4::ranef(fit)$GMA[,1], y=truth$GMAs, las=1,
##'        xlab="BLUP(GMA)", ylab="GMA")
##'   abline(a=0, b=1)
##' }
##' }
##' @noRd
lmerZ <- function(formula, data = NULL, REML = TRUE, control = lme4::lmerControl(),
                  start = NULL, verbose = 0L, subset, weights, na.action, offset,
                  contrasts = NULL, devFunOnly = FALSE, myblist) {
  stopifnot(
    is.list(myblist),
    !is.null(names(myblist)),
    all(names(myblist) %in% colnames(data)),
    all(sapply(myblist, function(x) {
      all(names(x) == c("ff", "sm", "nl", "cnms"))
    })),
    all(sapply(myblist, function(x) {
      is.integer(x$nl)
    }))
  )
  mc <- mcout <- match.call()
  missCtrl <- missing(control)
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if (!is.list(control)) {
      stop("'control' is not a list; use lmerControl()")
    }
    warning("passing control as list is deprecated: please use lmerControl() instead",
      immediate. = TRUE
    )
    control <- do.call(lme4::lmerControl, control)
  }
  mc$control <- control
  ## mc[[1]] <- quote(myLFormula)
  mc[[1]] <- myLFormula
  lmod <- eval(mc, parent.frame(1L))
  mcout$formula <- lmod$formula
  lmod$formula <- NULL
  devfun <- do.call(lme4::mkLmerDevfun, c(lmod, list(
    start = start,
    verbose = verbose, control = control
  )))
  if (devFunOnly) {
    return(devfun)
  }
  if (identical(control$optimizer, "none")) {
    stop("deprecated use of optimizer=='none'; use NULL instead")
  }
  opt <- if (length(control$optimizer) == 0) {
    s <- getStart(start, environment(devfun)$pp)
    list(par = s, fval = devfun(s), conv = 1000, message = "no optimization")
  } else {
    lme4::optimizeLmer(devfun,
      optimizer = control$optimizer, restart_edge = control$restart_edge,
      boundary.tol = control$boundary.tol, control = control$optCtrl,
      verbose = verbose, start = start, calc.derivs = control$calc.derivs,
      use.last.params = control$use.last.params
    )
  }
  cc <- lme4::checkConv(attr(opt, "derivs"), opt$par,
    ctrl = control$checkConv,
    lbound = environment(devfun)$lower
  )
  lme4::mkMerMod(environment(devfun), opt, lmod$reTrms,
    fr = lmod$fr,
    mc = mcout, lme4conv = cc
  )
}


myLFormula <- function(
  formula, data = NULL, REML = TRUE, subset, weights,
  na.action, offset, contrasts = NULL, control = lme4::lmerControl(), myblist,
  ...
) {
  control <- control$checkControl
  mf <- mc <- match.call()
  dontChk <- c("start", "verbose", "devFunOnly")
  dots <- list(...)
  do.call(checkArgs, c(list("lmer"), dots[!names(dots) %in%
    dontChk]))
  if (!is.null(dots[["family"]])) {
    mc[[1]] <- quote(lme4::glFormula)
    if (missing(control)) {
      mc[["control"]] <- lme4::glmerControl()
    }
    return(eval(mc, parent.frame()))
  }
  cstr <- "check.formula.LHS"
  checkCtrlLevels(cstr, control[[cstr]])
  denv <- checkFormulaData(formula, data, checkLHS = control$check.formula.LHS ==
    "stop")
  formula <- stats::as.formula(formula, env = denv)
  RHSForm(formula) <- lme4::expandDoubleVerts(RHSForm(formula))
  mc$formula <- formula
  m <- match(
    c("data", "subset", "weights", "na.action", "offset"),
    names(mf), 0L
  )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  fr.form <- lme4::subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) {
      assign(i, get(i, parent.frame()), environment(fr.form))
    }
  }
  mf$formula <- fr.form
  fr <- eval(mf, parent.frame())
  if (nrow(fr) == 0L) {
    stop("0 (non-NA) cases")
  }
  fr <- lme4::factorize(fr.form, fr, char.only = TRUE)
  attr(fr, "formula") <- formula
  attr(fr, "offset") <- mf$offset
  n <- nrow(fr)
  reTrms <- myMkReTrms(lme4::findbars(RHSForm(formula)), fr,
    myblist = myblist
  )
  wmsgNlev <- checkNlevels(reTrms$flist, n = n, control)
  wmsgZdims <- checkZdims(reTrms$Ztlist, n = n, control, allow.n = FALSE)
  if (anyNA(reTrms$Zt)) {
    stop(
      "NA in Z (random-effects model matrix): ", "please use ",
      shQuote("na.action='na.omit'"), " or ", shQuote("na.action='na.exclude'")
    )
  }
  wmsgZrank <- checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06)
  fixedform <- formula
  RHSForm(fixedform) <- lme4::nobars(RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(
    fixedfr,
    "terms"
  ), "predvars")
  attr(attr(fr, "terms"), "varnames.fixed") <- names(fixedfr)
  ranform <- formula
  RHSForm(ranform) <- lme4::subbars(RHSForm(reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(
    stats::terms(ranfr),
    "predvars"
  )
  X <- stats::model.matrix(fixedform, fr, contrasts)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lme4::lmerControl)[["check.rankX"]])[[1]]
  }
  X <- chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lme4::lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- checkScaleX(X, kind = scaleX.chk)
  list(
    fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula,
    wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)
  )
}


myMkReTrms <- function(
  bars, fr, myblist, drop.unused.levels = TRUE, reorder.terms = TRUE,
  reorder.vars = FALSE
) {
  if (!length(bars)) {
    stop("No random effects terms specified in formula",
      call. = FALSE
    )
  }
  stopifnot(is.list(bars), vapply(bars, is.language, NA), inherits(
    fr,
    "data.frame"
  ))
  names(bars) <- barnames(bars)
  term.names <- vapply(bars, deparse1, "")
  blist <- lapply(bars, mkBlist, fr, drop.unused.levels, reorder.vars = reorder.vars) # one component per random term
  ## begining of the code added to the official lme4 code ------------
  stopifnot(
    is.list(myblist),
    length(blist) == length(myblist),
    all(names(blist) == names(myblist))
  )
  for (i in seq_along(blist)) {
    if (is.na(myblist[[i]]$nl)) { # if NA, keep the element of blist as is
      next
    }
    stopifnot(all(names(myblist[[i]]) == c("ff", "sm", "nl", "cnms")))
    blist[[i]] <- myblist[[i]]
  }
  ## end of the code added to the official lme4 code -----------------
  nl <- vapply(blist, `[[`, 0L, "nl")
  if (reorder.terms) {
    if (any(diff(nl) > 0)) {
      ord <- rev(order(nl))
      blist <- blist[ord]
      nl <- nl[ord]
      term.names <- term.names[ord]
    }
  }
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(rbind, Ztlist)
  names(Ztlist) <- term.names
  q <- nrow(Zt)
  cnms <- lapply(blist, `[[`, "cnms")
  nc <- lengths(cnms)
  nth <- as.integer((nc * (nc + 1)) / 2)
  nb <- nc * nl
  if (sum(nb) != q) {
    stop(sprintf(
      "total number of RE (%d) not equal to nrow(Zt) (%d)",
      sum(nb), q
    ))
  }
  boff <- cumsum(c(0L, nb))
  thoff <- cumsum(c(0L, nth))
  Lambdat <- Matrix::t(do.call(Matrix::sparseMatrix, do.call(rbind, lapply(
    seq_along(blist),
    function(i) {
      mm <- matrix(seq_len(nb[i]), ncol = nc[i], byrow = TRUE)
      dd <- diag(nc[i])
      ltri <- lower.tri(dd, diag = TRUE)
      ii <- row(dd)[ltri]
      jj <- col(dd)[ltri]
      data.frame(i = as.vector(mm[, ii]) + boff[i], j = as.vector(mm[
        ,
        jj
      ]) + boff[i], x = as.double(rep.int(
        seq_along(ii),
        rep.int(nl[i], length(ii))
      ) + thoff[i]))
    }
  ))))
  thet <- numeric(sum(nth))
  ll <- list(
    Zt = Matrix::drop0(Zt), theta = thet, Lind = as.integer(Lambdat@x),
    Gp = unname(c(0L, cumsum(nb)))
  )
  ll$lower <- -Inf * (thet + 1)
  ll$lower[unique(Matrix::diag(Lambdat))] <- 0
  ll$theta[] <- is.finite(ll$lower)
  Lambdat@x[] <- ll$theta[ll$Lind]
  ll$Lambdat <- Lambdat
  fl <- lapply(blist, `[[`, "ff")
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else {
    asgn <- seq_along(fl)
  }
  names(fl) <- ufn
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist
  ll
}
