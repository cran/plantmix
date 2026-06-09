library(plantmix)
context("utils_stat")

test_that("normBiasError_percFALSE", {
  estim_vals <- c(0.1, 1 - 0.1, 2 + 0.3, 5, 6, 8, 9)
  true_vals <- c(0, 1, 2, 5, 6, 8, 9)
  expected <- (estim_vals - true_vals) / true_vals
  observed <- normBiasError(estim_vals, true_vals, perc = FALSE)
  expect_equal(observed, expected)
})

test_that("normBiasError_percTRUE", {
  estim_vals <- c(0.1, 1 - 0.1, 2 + 0.3, 5, 6, 8, 9)
  true_vals <- c(0, 1, 2, 5, 6, 8, 9)
  expected <- 100 * (estim_vals - true_vals) / true_vals
  observed <- normBiasError(estim_vals, true_vals)
  expect_equal(observed, expected)
})

test_that("jaccard_numeric", {
  s1 <- c(0, 1, 2, 5, 6, 8, 9)
  s2 <- c(0, 2, 3, 4, 5, 7, 9)
  card_inter <- 4
  card_union <- 10
  expected <- card_inter / card_union
  observed <- jaccard(s1, s2)
  expect_equal(observed, expected)
})

test_that("jaccard_character", {
  s1 <- c("geno1", "geno2")
  s2 <- c("geno1", "geno3")
  card_inter <- 1
  card_union <- 3
  expected <- card_inter / card_union
  observed <- jaccard(s1, s2)
  expect_equal(observed, expected)
})

test_that("jaccardBtwCombs", {
  combs <- list(
    mix1 = c("comp1", "comp2", "comp3"),
    mix2 = c("comp1", "comp3"),
    mix3 = c("comp2", "comp3")
  )
  j12 <- jaccard(combs$mix1, combs$mix2)
  j13 <- jaccard(combs$mix1, combs$mix3)
  j23 <- jaccard(combs$mix2, combs$mix3)
  expected <- matrix(
    c(
      1, j12, j13,
      j12, 1, j23,
      j13, j23, 1
    ),
    nrow = 3, ncol = 3, byrow = TRUE,
    dimnames = list(names(combs), names(combs))
  )
  observed <- jaccardBtwCombs(combs)
  expect_is(observed, "sparseMatrix")
  observed <- as.matrix(observed)
  expect_equal(observed, expected)
})

test_that("getShannonIndex", {
  p <- c(0.2, 0.8)
  expected <- -(p[1] * log(p[1]) + p[2] * log(p[2]))
  observed <- getShannonIndex(p)
  expect_equal(observed, expected)
})

test_that("getGenoPropsFromCombs", {
  combs <- list(
    "comb1" = c("comp1", "comp2"),
    "comb2" = c("comp1", "comp3")
  )
  expected <- c(
    "comp1" = 2 / 4,
    "comp2" = 1 / 4,
    "comp3" = 1 / 4
  )
  observed <- getGenoPropsFromCombs(combs)
  expect_equal(observed, expected)
})
