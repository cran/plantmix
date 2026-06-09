##' Matrix-variate Normal distribution
##'
##' Random generation for the matrix-variate Normal distribution.
##' See \url{https://en.wikipedia.org/wiki/Matrix_normal_distribution}.
##' @param n number of observations
##' @param M mean matrix; if missing, will be replaced by a matrix of zeroes; if \code{dimnames(M)} is not NULL, the column and row names will be compared to those of U and V (if any) and, in the absence of conflict, propagated to the output
##' @param U between-row covariance matrix
##' @param V between-column covariance matrix
##' @param pivot 2-element vector with values TRUE/FALSE/"auto", where TRUE (FALSE) means using pivoting (or not) for Choleski decomposition of U and/or V (see \code{\link[base]{chol}}); useful when U and/or V are singular; with "auto", this will be automatically determined
##' @return array
##' @author Timothee Flutre
##' @examples
##' set.seed(1859)
##' Sigma <- matrix(c(3,2,2,4), nrow=2, ncol=2)
##' rho <- Sigma[2,1] / prod(sqrt(diag(Sigma)))
##' samples <- rmatnorm(n=100, M=matrix(0, nrow=10^3, ncol=2),
##'                     U=diag(10^3), V=Sigma)
##' tmp <- t(apply(samples, 3, function(mat){
##'   c(var(mat[,1]), var(mat[,2]), cor(mat[,1], mat[,2]))
##' }))
##' summary(tmp) # corresponds well to Sigma
##' @export
rmatnorm <- function(n = 1, M, U, V, pivot = c(U = "auto", V = "auto")) {
  if (missing(M)) {
    M <- matrix(0, nrow = nrow(U), ncol = ncol(V))
    dimnames(M) <- list(rownames(U), colnames(V))
  }
  stopifnot(
    is.matrix(M),
    is.matrix(U),
    is.matrix(V),
    nrow(M) == nrow(U),
    ncol(M) == nrow(V),
    nrow(U) == ncol(U),
    nrow(V) == ncol(V),
    is.vector(pivot),
    all(c("U", "V") %in% names(pivot)),
    all(pivot %in% c(TRUE, FALSE, "auto"))
  )
  if (!is.null(colnames(M))) {
    if (!is.null(colnames(V))) {
      stopifnot(all(colnames(M) == colnames(V)))
    }
    if (!is.null(rownames(V))) {
      stopifnot(all(colnames(M) == rownames(V)))
    }
  }
  if (!is.null(rownames(M))) {
    if (!is.null(colnames(U))) {
      stopifnot(all(rownames(M) == colnames(U)))
    }
    if (!is.null(rownames(U))) {
      stopifnot(all(rownames(M) == rownames(U)))
    }
  }

  ## chol() returns upper triangular factor of Cholesky decomp
  if (pivot["U"] == "auto") {
    chol.tU <- chol(t(U), pivot = isSingular(U)) # A
  } else {
    chol.tU <- chol(t(U), pivot = pivot["U"])
  } # A
  if (pivot["V"] == "auto") {
    chol.V <- chol(V, pivot = isSingular(V)) # B
  } else {
    chol.V <- chol(V, pivot = pivot["V"])
  } # B

  ## for X ~ MN(M, AA', B'B): draw Z ~ MN(0, I, I), then X = M + A Z B
  tmp <- lapply(1:n, function(i) {
    Z <- matrix(
      data = stats::rnorm(n = nrow(M) * ncol(M), mean = 0, sd = 1),
      nrow = nrow(M), ncol = ncol(M)
    )
    matrix(
      data = M + chol.tU %*% Z %*% chol.V,
      nrow = nrow(M), ncol = ncol(M)
    )
  })

  out <- array(
    data = do.call(c, tmp),
    dim = c(nrow(M), ncol(M), n)
  )
  dimnames(out) <- list(rownames(M), colnames(M), NULL)

  return(out)
}
