## from https://github.com/timflutre/rutilstimflutre/blob/master/R/quantgen.R

##' Simulate SNP genotypes
##'
##' Simulates SNP genotypes as allele dose additively encoded, i.e. 0,1,2, using correlated allele frequencies to mimick genetic structure.
##' @param nb_genos number of genotypes (i.e. individuals) per population
##' @param nb_snps number of SNPs
##' @param div_pops matrix of divergence among populations, with a diagonal of 1's; the closer off-diagonal values are from 1; the weaker the divergence; the further, the stronger
##' @param geno_IDs vector of genotype identifiers (if NULL, will be "geno001", etc)
##' @param snp_IDs vector of SNP identifiers (if NULL, will be "snp001", etc)
##' @return matrix with genotypes in rows and SNPs in columns
##' @seealso \code{\link{estimGRM}}
##' @author Timothee Flutre thanks to code from Andres Legarra
##' @examples
##' ## weak divergences among populations:
##' weak_div_pops <- diag(3)
##' weak_div_pops[upper.tri(weak_div_pops)] <- 0.9
##' weak_div_pops[lower.tri(weak_div_pops)] <- weak_div_pops[upper.tri(weak_div_pops)]
##' weak_div_pops
##'
##' ## strong divergences among populations:
##' strong_div_pops <- diag(3)
##' strong_div_pops[upper.tri(strong_div_pops)] <- 0.5
##' strong_div_pops[lower.tri(strong_div_pops)] <- strong_div_pops[upper.tri(strong_div_pops)]
##' strong_div_pops
##'
##' M <- simulGenosDoseStruct(div_pops=weak_div_pops)
##' A <- estimGRM(M)
##' imageMat(A, "Weak divergence")
##'
##' M <- simulGenosDoseStruct(div_pops=strong_div_pops)
##' A <- estimGRM(M)
##' imageMat(A, "Strong divergence")
##' @export
simulGenosDoseStruct <- function(nb_genos = c(100, 120, 80),
                                 nb_snps = 1000,
                                 div_pops = diag(3) * 0.5 + 0.5,
                                 geno_IDs = NULL, snp_IDs = NULL) {
  stopifnot(
    requireNamespace("MASS", quietly = TRUE),
    length(nb_genos) > 1,
    all(nb_genos %% 1 == 0),
    length(nb_genos) == nrow(div_pops),
    ncol(div_pops) == nrow(div_pops),
    all(sapply(diag(div_pops), all.equal, 1)),
    all(c(div_pops) >= 0),
    all(c(div_pops) <= 1)
  )
  if (!is.null(geno_IDs)) {
    stopifnot(length(geno_IDs) == sum(nb_genos))
  }
  if (!is.null(snp_IDs)) {
    stopifnot(length(snp_IDs) == nb_snps)
  }

  nbPops <- length(nb_genos)
  idx.pops <- matrix(nrow = nbPops, ncol = 2)
  idx.pops[1, ] <- c(1, nb_genos[1])
  for (j in 2:nbPops) {
    idx.pops[j, ] <- c(
      sum(nb_genos[1:(j - 1)]) + 1,
      sum(nb_genos[1:(j - 1)]) + nb_genos[j]
    )
  }

  ## draw allele frequencies from "not-too-different populations" with correlated allele frequencies
  allFreqs <- MASS::mvrnorm(n = nb_snps, mu = rep(0, nbPops), Sigma = div_pops)
  allFreqs <- qbeta(pnorm(allFreqs), 2, 2)

  ## sample SNP genotypes and assign {0,1,2} coding
  M <- lapply(1:nb_snps, function(i) {
    tmp <- rep(NA, sum(nb_genos))
    for (j in 1:nbPops) {
      idx <- idx.pops[j, 1]:idx.pops[j, 2]
      tmp[idx] <- sample(0:1, nb_genos[j],
        replace = TRUE,
        prob = c(1 - allFreqs[i, j], allFreqs[i, j])
      ) +
        sample(0:1, nb_genos[j],
          replace = TRUE,
          prob = c(1 - allFreqs[i, j], allFreqs[i, j])
        )
    }
    tmp
  })
  M <- do.call(cbind, M)

  if (is.null(geno_IDs)) {
    geno_IDs <- sprintf(
      fmt = paste0("geno%0", floor(log10(nb_genos)) + 1, "i"),
      1:sum(nb_genos)
    )
  }
  rownames(M) <- geno_IDs

  if (is.null(snp_IDs)) {
    snp_IDs <- sprintf(
      fmt = paste0("snp%0", floor(log10(nb_snps)) + 1, "i"),
      1:nb_snps
    )
  }
  colnames(M) <- snp_IDs

  return(M)
}

##' Estimate a GRM
##'
##' Estimates a genomic relationship matrix with the first estimator from VanRaden (2008).
##' @param M matrix of SNP genotypes, with individuals in rows and SNPs in columns, encoded in allele dose in \[0,2\]
##' @param AFs vector of allele frequencies; if NULL, will be estimated from \code{M}
##' @return matrix
##' @seealso \code{\link{simulGenosDoseStruct}}
##' @author Timothee Flutre
##' @examples
##' strong_div_pops <- diag(3)
##' strong_div_pops[upper.tri(strong_div_pops)] <- 0.5
##' strong_div_pops[lower.tri(strong_div_pops)] <- strong_div_pops[upper.tri(strong_div_pops)]
##' strong_div_pops
##' M <- simulGenosDoseStruct(div_pops=strong_div_pops)
##' A <- estimGRM(M)
##' imageMat(A, "Strong divergence")
##' @export
estimGRM <- function(M, AFs = NULL) {
  stopifnot(
    is.matrix(M),
    all(!is.na(M)),
    all(M >= 0),
    all(M <= 2)
  )
  if (!is.null(AFs)) {
    stopifnot(
      is.vector(AFs),
      is.numeric(AFs),
      all(!is.na(AFs)),
      all(AFs >= 0.0),
      all(AFs <= 1.0),
      length(AFs) == ncol(M)
    )
  }

  N <- nrow(M) # nb of individuals
  P <- ncol(M) # nb of SNPs

  if (is.null(AFs)) {
    AFs <- colMeans(M) / 2
  }

  ## implementation as in VanRaden (2008)
  ## M <- M - 1 # recode genotypes as {-1,0,1}
  ## Pmat <- matrix(rep(1, N)) %*% (2 * (afs - 0.5))
  ## Z <- M - Pmat

  ## implementation as in Vitezica et al (2013)
  tmp <- matrix(rep(1, N)) %*% (2 * AFs)
  Z <- M - tmp

  GRM <- tcrossprod(Z, Z) / (2 * sum(AFs * (1 - AFs)))

  return(GRM)
}
