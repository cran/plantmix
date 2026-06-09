##' Normalized bias error
##'
##' Returns the normalized bias error (NBE).
##' @param estim numeric vector of estimates
##' @param truth numeric vector of true values
##' @param perc logical; if TRUE, the return value will be in percentage
##' @return numeric vector
##' @author Timothee Flutre
##' @export
normBiasError <- function(estim, truth, perc = TRUE) {
  out <- (estim - truth) / truth
  if (perc) {
    100 * out
  } else {
    out
  }
}

##' Information criterion
##'
##' Returns an information criterion, among AIC, AICc and BIC.
##' Caution: choosing \code{k} and \code{n} may not be straightforward for certain models, such as mixed models.
##' @param k number of parameters, sometimes also called "effective degrees of freedom"
##' @param lnLmax value of the log-likelihood when maximized
##' @param n sample size
##' @param type specific information criterion to be returned
##' @return numeric with attributes \code{k} and \code{n} (if not NULL)
##' @author Timothee Flutre
##' @export
infoCriterion <- function(k, lnLmax, n = NULL, type = "AIC") {
  stopifnot(
    is.numeric(k),
    is.numeric(lnLmax),
    type %in% c("AIC", "AICc", "BIC")
  )
  if (type %in% c("AICc", "BIC")) {
    stopifnot(
      !is.null(n),
      is.numeric(n)
    )
  }

  out <- -2 * lnLmax

  if (grepl("^AIC", type)) {
    out <- out + 2 * k
    if (type == "AICc") {
      out <- out + (2 * k^2 + 2 * k) / (n - k - 1)
    }
  } else if (type == "BIC") {
    out <- out + log(n) * k
  }

  attr(out, "k") <- k
  if (!is.null(n)) {
    attr(out, "n") <- n
  }

  return(out)
}

##' Jaccard index
##'
##' Returns the Jaccard index between two sets.
##' @examples
##' jaccard(c(0, 1, 2, 5, 6, 8, 9), c(0, 2, 3, 4, 5, 7, 9))
##' jaccard(c("a","b","c"), c("a","b"))
##' jaccard(c(3, 4), c(1,8))
##' @noRd
jaccard <- function(s1, s2) {
  out <- NA
  inter <- length(intersect(s1, s2))
  uni <- length(s1) + length(s2) - inter
  if (isTRUE(all.equal(inter / uni, 0.0))) {
    out <- 0.0
  } else {
    out <- inter / uni
  }
  return(out)
}

##' Returns a n x n sparse matrix from a list of n combinations (mixtures), the distance between two combinations being their Jaccard index.
##' @examples
##' combs <- list(mix1=c("comp1","comp2"),
##'               mix2=c("comp1","comp3"),
##'               mix3=c("comp2","comp3"))
##' (Kmix <- jaccardBtwCombs(combs))
##' combs <- list(mix1=c("comp1","comp2","comp3"),
##'               mix2=c("comp1","comp3"),
##'               mix3=c("comp2","comp3"))
##' (Kmix <- jaccardBtwCombs(combs))
##' @noRd
jaccardBtwCombs <- function(combs) {
  stopifnot(is.list(combs))

  if (FALSE) { # slow: two nested for loops
    out <- matrix(0,
      nrow = length(combs), ncol = length(combs),
      dimnames = list(names(combs), names(combs))
    )
    for (i in 1:(nrow(out) - 1)) {
      for (j in (i + 1):nrow(out)) {
        out[i, j] <- jaccard(combs[[i]], combs[[j]])
      }
    }
    diag(out) <- 1
    out <- symmetrize(out, to_copy = "upper")
  } else { # fast: vectorized code
    comps <- sort(unique(unlist(combs)))
    comp_index <- setNames(seq_len(length(comps)), comps)
    i <- rep(seq_len(length(combs)), lengths(combs))
    j <- unlist(lapply(combs, function(x) comp_index[x]))
    incid <- sparseMatrix( # incidence matrix
      i = i, j = j, x = 1, dims = c(length(combs), length(comps)),
      dimnames = list(names(combs), comps)
    )
    sizes <- rowSums(incid)
    inter_mat <- incid %*% t(incid) # intersections
    tmp <- summary(inter_mat)
    inter <- tmp$x
    uni <- sizes[tmp$i] + sizes[tmp$j] - tmp$x
    jacc_vals <- inter / uni
    out <- sparseMatrix(
      i = tmp$i,
      j = tmp$j,
      x = jacc_vals,
      dims = dim(inter_mat)
    )
  }
  dimnames(out) <- list(names(combs), names(combs))

  return(out)
}

##' @noRd
getShannonIndex <- function(p) {
  stopifnot(isTRUE(all.equal(sum(p), 1.0)))
  -sum(p * log(p))
}

##' @noRd
getGenoPropsFromCombs <- function(combs) {
  stopifnot(
    is.list(combs),
    all(sapply(combs, is.vector))
  )
  x <- table(unlist(combs))
  out <- x / sum(x)
  out <- setNames(as.vector(out), names(x))
  return(out)
}

##' @noRd
getOptimalProps <- function(levGenos, nbMixes) {
  nbGenos <- length(levGenos)
  out <- rep(NA, nbGenos)
  nbMixesPerGeno <- nbMixes / nbGenos
  if (nbMixesPerGeno %% 1 == 0) {
    out <- rep(nbMixesPerGeno, nbGenos)
  } else {
    out <- rep(floor(nbMixes / nbGenos), nbGenos)
    remainMixes <- nbMixes - sum(out)
    if (remainMixes <= length(levGenos)) {
      out[1:remainMixes] <- out[1:remainMixes] + 1
    } else {
      out[1] <- out[1] + nbMixes - sum(out)
    } # TODO: improve this
  }
  stopifnot(sum(out) == nbMixes)
  names(out) <- levGenos
  out <- out / sum(out)
  return(out)
}
