## dgTMatrix: class of sparse matrices in triplet form
sparseMat4Tmb <- function(mat) {
  as(as(as(mat, "dMatrix"), "generalMatrix"), "TsparseMatrix")
}

## @param capture output from \code{capture.output(nlminb(..., control=list("trace"=1)))}
traceFromNlminb <- function(capture) {
  stopifnot(is(capture, "character"))
  out <- NULL
  idx <- grep("^[ 0-9]+:", capture)
  if (length(idx) > 0) {
    out <- capture[idx]
    out <- strcapture(
      "^(\\s+[0-9]+):\\s+([.0-9]+):",
      out,
      data.frame(
        iter = numeric(),
        objfn = numeric()
      )
    )
  }
  return(out)
}

##' @noRd
prepParamBoot <- function(fit) {
  stopifnot(all(c("obj", "inputs4TMB") %in% names(fit)))

  out <- list()

  vec_estim_params <- fit$obj$env$last.par.best
  out$vec <- vec_estim_params

  list_estim_params <- list()
  for (param in names(fit$inputs4TMB$listParams)) {
    idx <- grep(paste0("^", param, "$"), names(vec_estim_params))

    if (length(idx) == 0) {
      stopifnot(param %in% names(fit$inputs4TMB$map))
      list_estim_params[[param]] <- fit$inputs4TMB$listParams[[param]]
    } else {
      if (is.vector(fit$inputs4TMB$listParams[[param]])) {
        list_estim_params[[param]] <- vec_estim_params[idx]
        names(list_estim_params[[param]]) <- NULL
      } else if (is.matrix(fit$inputs4TMB$listParams[[param]])) {
        list_estim_params[[param]] <-
          matrix(vec_estim_params[idx],
            ncol = ncol(fit$inputs4TMB$listParams[[param]]),
            byrow = FALSE
          )
      } else {
        msg <- paste0(param, "is neither a vector nor a matrix")
        stop(msg)
      }
    }
  }
  out$list <- list_estim_params

  return(out)
}

##' Parametric bootstrap with TMB
##'
##' Performs parametric bootstrap with TMB and \code{\link{nlminb}} on a fitted model to allow uncertainty quantification.
##' @param fit list corresponding to a model fitted with TMB
##' @param nb_boot number of parametric bootstraps to perform
##' @param mc.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously (see the "parallel" package)
##' @return list with as many components as \code{nb_boot}
##' @author Jemay Salomon, Timothee Flutre
##' @examples
##' ## see the example of `fitDBVSBVinter`
##' ## then run (slow if nb_boot is high!): paramBoot4TMB(fitTmb)
##' @export
paramBoot4TMB <- function(fit, nb_boot = 5, mc.cores = 1) {
  stopifnot(all(c("obj", "inputs4TMB") %in% names(fit)))
  if (mc.cores > 1) {
    msg <- paste0("not yet parallelized, the value of mc.cores is ignored")
    warning(msg)
  }

  estim_params <- prepParamBoot(fit)

  ## TODO: parallelize with parallel::parLapply()
  ## https://stackoverflow.com/a/19281611/597069
  ## Inspired by https://kaskr.github.io/adcomp/Simulation.html
  out <- replicate(nb_boot,
    {
      data_sim <- fit$obj$simulate(par = estim_params$vec, complete = TRUE)
      resp_var_names_sim <- grep("_sim$", names(data_sim), value = TRUE)
      for (rvns in resp_var_names_sim) {
        resp_var_name <- sub("_sim$", "", rvns)
        data_sim[[resp_var_name]] <- data_sim[[rvns]]
        data_sim[[rvns]] <- NULL
      }
      f_boot <- MakeADFun(
        data = data_sim,
        parameters = estim_params$list,
        map = fit$inputs4TMB$map,
        random = fit$inputs4TMB$vecRnd,
        DLL = "plantmix_TMBExports",
        silent = TRUE
      )
      nlminb(f_boot$par, f_boot$fn, f_boot$gr)$par
    },
    simplify = FALSE
  )

  return(out)
}
