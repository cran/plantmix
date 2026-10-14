library(plantmix)
context("uncertainty_RYM")

test_that("RYMperRep_weighted", {
  dat <- data.frame(
    ID = rep(c("geno1", "geno2", "mixg1g2"), times = 2),
    block = c("A", "A", "A", "B", "B", "B"),
    yield = c(51, 45, 39, 43, 25, 50),
    stringsAsFactors = TRUE
  )
  expected <- list(
    "A" = data.frame(
      RYM = c(1, 1, 39 / c(0.6 * 51 + 0.4 * 45)),
      block = c("A", "A", "A"),
      ID = c("geno1", "geno2", "mixg1g2"),
      stringsAsFactors = TRUE
    ),
    "B" = data.frame(
      RYM = c(1, 1, 50 / c(0.6 * 43 + 0.4 * 25)),
      block = c("B", "B", "B"),
      ID = c("geno1", "geno2", "mixg1g2"),
      stringsAsFactors = TRUE
    )
  )
  expected <- do.call(rbind, expected)
  rownames(expected) <- NULL
  observed <- RYMperRep(
    dat, "yield", "ID", "block",
    list(
      "geno1" = c("geno1" = 1.0),
      "geno2" = c("geno2" = 1.0),
      "mixg1g2" = c("geno1" = 0.6, "geno2" = 0.4)
    )
  )
  expect_equal(observed, expected)
})
