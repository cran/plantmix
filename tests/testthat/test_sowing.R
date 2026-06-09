library(plantmix)
context("sowing")

test_that("nbSeedsToSownInPure", {
  sowingArea <- 8.4
  sowingDensity <- 160

  expected <- 1.344

  observed <- nbSeedsToSownInPure(sowingArea, sowingDensity)

  expect_equal(observed, expected)
})

test_that("mixSowingWeight", {
  (tmp <- nbSeedsToSownInPure(sowingArea = 9.52, sowingDensity = 160))
  TKWs <- c(
    "var1" = 38.605,
    "var2" = 40.051,
    "var6" = 36.251,
    "var8" = 33.368
  )
  (pureSowingWeights <- TKWs * tmp)
  stands <- c(
    "var1" = "var1",
    "mix8" = "var1-var2",
    "mix34" = "var1-var8",
    "mix50" = "var1-var6-var8"
  )

  expected <- list(
    "var1" = pureSowingWeights["var1"],
    "mix8" = c("var1" = 29.401568, "var2" = 30.502842),
    "mix34" = c("var1" = 29.401568, "var8" = 25.413069),
    "mix50" = c("var1" = 19.405035, "var6" = 18.221783, "var8" = 16.772625)
  )

  observed <- mixSowingWeight(pureSowingWeights, stands)

  expect_equal(observed, expected, tolerance = 0.01)
})
