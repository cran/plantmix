library(plantmix)
context("uncertainty_RY")

test_that("RYperRep", {
  dat <- data.frame(
    ID = c(
      "geno1", "geno1",
      "geno2", "geno2",
      "mixg1g2", "mixg1g2", "mixg1g2", "mixg1g2"
    ),
    focal = c(
      "geno1", "geno1",
      "geno2", "geno2",
      "geno1", "geno2",
      "geno1", "geno2"
    ),
    prop = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
    block = c("A", "B", "A", "B", "A", "A", "B", "B"),
    yield = c(51, 45, 39, 43, 25, 22, 26, 25)
  )
  expected <- data.frame(
    ID = rep("mixg1g2", 4),
    focal = c("geno1", "geno2", "geno1", "geno2"),
    block = c("A", "A", "B", "B"),
    RY = c(25 / 51, 22 / 39, 26 / 45, 25 / 43),
    stringsAsFactors = TRUE
  )
  observed <- RYperRep(
    dat, "yield", "ID", "focal", "block",
    list("mixg1g2" = c("geno1", "geno2"))
  )
  expect_equal(observed, expected)
})
