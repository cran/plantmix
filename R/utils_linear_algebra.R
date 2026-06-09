##' Image of a matrix
##'
##' Plots an image of a matrix.
##' @param mat matrix
##' @param title optional title of the plot; if missing, it will be the matrix name and dimensions
##' @param col optional vector of colors; if missing, it will be a grey color gradient
##' @param breaks optional vector of breaks to go with col
##' @return nothing
##' @author Timothee Flutre
##' @examples
##' M <- diag(c(rep(2, 3), rep(15, 5)))
##' imageMat(M) # see how the the "2" values are ignored!
##' imageMat(M, col=c("lightgrey","red","blue"), breaks=c(-1,1,3,20))
##' @export
imageMat <- function(mat, title, col, breaks) {
  if (is(mat, "Matrix")) {
    msg <- paste0(
      "your matrix is in the Matrix format",
      "; you should rather use the 'image' function of this package"
    )
    warning(msg, immediate. = TRUE)
    mat <- as.matrix(mat)
  }
  if (missing(title)) {
    title <- deparse(substitute(mat))
    title <- paste0(title, " (", nrow(mat), " x ", ncol(mat), ")")
  }
  if (missing(col)) {
    col <- grey.colors(length(table(c(mat))), rev = TRUE)
  }
  image(t(mat)[, nrow(mat):1], axes = FALSE, main = title, col = col, breaks = breaks)
}

##' @noRd
symmetrize <- function(mat, to_copy = "upper") {
  stopifnot(
    is.matrix(mat),
    nrow(mat) == ncol(mat),
    to_copy %in% c("upper", "lower")
  )
  out <- mat
  if (to_copy == "upper") {
    out[lower.tri(out)] <- t(out)[lower.tri(out)]
  } else {
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
  }
  stopifnot(isSymmetric(out))
  return(out)
}

##' Singular matrix
##'
##' Assesses if a matrix is singular, i.e. not full rank, by comparing its condition number with the relative machine precision.
##' As \href{http://www.stat.wisc.edu/~st849-1/Rnotes/ModelMatrices.html#sec-2_2}{explained by Douglas Bates}, "a matrix is regarded as being numerically singular when its reciprocal condition number is less than the relative machine precision".
##' @param x matrix
##' @return logical
##' @author Timothee Flutre
##' @noRd
isSingular <- function(x) {
  return(1 / kappa(x) <= 100 * .Machine$double.eps)
}

##' Vec operator
##'
##' Applies the vec operator to a matrix, i.e., concatenate its columns, and save the names, too, if any.
##' @param mat matrix
##' @param sep separator of column and row names; used only if \code{mat} has dimnames
##' @return vector, possibly with names as <rowname><sep><colname>
##' @author Timothee Flutre
##' @seealso \code{\link{invvec}}
##' @examples
##' mat <- matrix(c(1, 2, 3,
##'                 4, 5, 6),
##'               nrow = 2, ncol = 3, byrow = TRUE,
##'               dimnames = list(as.character(1:2), letters[1:3]))
##' mat
##' vec(mat)
##' @export
vec <- function(mat, sep = "-") {
  out <- c(mat)
  if (!is.null(dimnames(mat))) {
    names(out) <- paste0(
      rownames(mat),
      sep,
      rep(colnames(mat), each = nrow(mat))
    )
  }
  return(out)
}

##' Inverse vec operator
##'
##' Applies the inverse vec operator to a vector.
##' @param x vector
##' @param n_row number of rows of the output matrix
##' @param n_col number of columns of the output matrix
##' @param sep separator used to retrieve row and column names, only if \code{x} has names
##' @return matrix
##' @seealso \code{\link{vec}}
##' @author Timothee Flutre
##' @examples
##' mat <- matrix(c(1, 2, 3,
##'                 4, 5, 6),
##'               nrow = 2, ncol = 3, byrow = TRUE,
##'               dimnames = list(as.character(1:2), letters[1:3]))
##' mat
##' (x <- vec(mat))
##' invvec(x, n_col = 3, sep = "-")
##' @export
invvec <- function(x, n_col, n_row = NULL, sep = NULL) {
  if (is.null(n_row)) {
    n_row <- length(x) / n_col
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
    stopifnot(is.wholenumber(n_row))
  }
  out <- matrix(x, nrow = n_row, ncol = n_col, byrow = FALSE)
  if (FALSE) {
    ## https://math.stackexchange.com/a/3122442
    I_n_row <- diag(n_row)
    I_n_col <- diag(n_col)
    out2 <- (t(vec(I_n_col)) %x% I_n_row) %*% (I_n_col %x% x)
  }
  if (all(!is.null(names(x)), !is.null(sep))) {
    tmp <- strsplit(names(x), paste0("\\", sep))
    rowNs <- unique(sapply(tmp, `[`, 1))
    colNs <- unique(sapply(tmp, `[`, 2))
    rownames(out) <- rowNs
    colnames(out) <- colNs
  }
  return(out)
}

##' Inverse vec for mixtures
##'
##' Tranforms a vector of component observations into a matrix with one row per mixture and components in colums.
##' Requires all mixtures to have the same number of components, and all the obervationss from the same mixture to be one after the others, in the same order for each mixture.
##' @param x vector; make sure that it is correctly sorted beforehand
##' @param nb_comps number of components in each mixture
##' @return matrix with \code{nb_comps} columns
##' @author Timothee Flutre
##' @examples
##' datL <- data.frame(ID = c("g1+g2", "g1+g2", "g1+g2", "g1+g2"),
##'   focal = c("g1", "g2", "g1", "g2"),
##'   neighbor = c("g2", "g1", "g2", "g1"),
##'   block = c("A", "A", "B", "B"),
##'   yield = c(1, 2, 3, 4)
##' )
##' invvecMixes(datL$yield, 2)
##' @export
invvecMixes <- function(x, nb_comps) {
  if (is(x, "matrix")) {
    stopifnot(ncol(x) == 1)
    x <- x[, 1]
  }
  stopifnot(
    is.vector(x),
    nb_comps >= 2
  )

  n <- length(x)
  n_row <- n / nb_comps
  out <- matrix(NA, nrow = n_row, ncol = nb_comps)
  I_row <- diag(n_row)
  for (j in 1:nb_comps) {
    idx <- seq(j, n, by = nb_comps)
    out[, j] <- I_row %*% x[idx]
  }
  return(out)
}
