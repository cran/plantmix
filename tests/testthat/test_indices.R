library(plantmix)
context("indices")

test_that("RYP", {
  ## example of a binary varietal mixture with unequal proportions:
  dat <- data.frame(
    ID = c("geno1", "geno2", "mixg1g2", "mixg1g2"),
    focal = c("geno1", "geno2", "geno1", "geno2"),
    prop = c(1, 1, 0.6, 0.4),
    yield = c(50, 40, 25, 22)
  )
  expected <- data.frame(dat,
    RYP = c(1, 1, 25 / (0.6 * 50), 22 / (0.4 * 40))
  )
  observed <- RYP(dat)
  expect_equal(observed, expected)
})

test_that("RY", {
  ## example of a binary varietal mixture with unequal proportions:
  dat <- data.frame(
    ID = c("geno1", "geno2", "mixg1g2", "mixg1g2"),
    focal = c("geno1", "geno2", "geno1", "geno2"),
    prop = c(1, 1, 0.6, 0.4),
    yield = c(50, 40, 25, 22)
  )
  expected <- data.frame(dat,
    RY = c(1, 1, 25 / 50, 22 / 40)
  )
  observed <- RY(dat)
  expect_equal(observed, expected)
})

test_that("RYT", {
  ## example of a binary varietal mixture with unequal proportions:
  dat <- data.frame(
    ID = c("geno1", "geno2", "mixg1g2", "mixg1g2"),
    focal = c("geno1", "geno2", "geno1", "geno2"),
    prop = c(1, 1, 0.6, 0.4),
    yield = c(50, 40, 25, 22)
  )
  expected <- data.frame(RY(dat),
    RYT = c(NA, NA, (25 / 50) + (22 / 40), (25 / 50) + (22 / 40))
  )
  observed <- RYT(dat)
  expect_equal(observed, expected)
})

test_that("RYM_withStandIDs_withProps", {
  ## example of a binary varietal mixture with unequal proportions:
  dat <- data.frame(
    ID = c("geno1", "geno2", "mixg1g2"),
    comps = c("geno1", "geno2", "geno1-geno2"),
    props = c("1", "1", "0.6-0.4"),
    yield = c(50, 40, 47)
  )
  expected <- data.frame(dat,
    RYM = c(NA, NA, 47 / (0.6 * 50 + 0.4 * 40))
  )
  observed <- RYM(dat,
    colIDstand = "ID", colIDcomps = "comps", colProps = "props",
    sep = "-"
  )
  expect_equal(observed, expected)
})

test_that("RYM_withStandIDs_withProps_plusSign", {
  ## example of a binary varietal mixture with unequal proportions:
  dat <- data.frame(
    ID = c("geno1", "geno2", "mixg1g2"),
    comps = c("geno1", "geno2", "geno1+geno2"),
    props = c("1", "1", "0.6+0.4"),
    yield = c(50, 40, 47)
  )
  expected <- data.frame(dat,
    RYM = c(NA, NA, 47 / (0.6 * 50 + 0.4 * 40))
  )
  observed <- RYM(dat,
    colIDstand = "ID", colIDcomps = "comps", colProps = "props",
    sep = "+"
  )
  expect_equal(observed, expected)
})

test_that("RYM_withoutStandIDs_withProps", {
  ## example of a binary varietal mixture with unequal proportions:
  dat <- data.frame(
    comps = c("geno1", "geno2", "geno1-geno2"),
    props = c("1", "1", "0.6-0.4"),
    yield = c(50, 40, 47)
  )
  expected <- data.frame(dat,
    RYM = c(NA, NA, 47 / (0.6 * 50 + 0.4 * 40))
  )
  observed <- RYM(dat, colIDcomps = "comps", colProps = "props", sep = "-")
  expect_equal(observed, expected)
})

test_that("RYM_withoutStandIDs_withProps_plusSign", {
  ## example of a binary varietal mixture with unequal proportions:
  dat <- data.frame(
    comps = c("geno1", "geno2", "geno1+geno2"),
    props = c("1", "1", "0.6+0.4"),
    yield = c(50, 40, 47)
  )
  expected <- data.frame(dat,
    RYM = c(NA, NA, 47 / (0.6 * 50 + 0.4 * 40))
  )
  observed <- RYM(dat, colIDcomps = "comps", colProps = "props", sep = "+")
  expect_equal(observed, expected)
})

test_that("RYM_withoutStandIDs_withoutProps", {
  ## example of a binary varietal mixture with unequal proportions:
  dat <- data.frame(
    comps = c("geno1", "geno2", "geno1-geno2"),
    yield = c(50, 40, 47)
  )
  expected <- data.frame(dat,
    RYM = c(NA, NA, 47 / (0.5 * 50 + 0.5 * 40))
  )
  observed <- RYM(dat, sep = "-")
  expect_equal(observed, expected)
})

test_that("RY_RD18", {
  perfMix <- setNames(
    c(40, 50),
    c("mix2", "mix1")
  )
  perfPur <- setNames(
    c(70, 35, 20),
    c("varA1", "varC2", "varB3")
  )
  mix2pur <- list(
    "mix1" = setNames(c(0.5, 0.5), c("varB3", "varC2")),
    "mix2" = setNames(c(0.5, 0.2, 0.3), c("varC2", "varA1", "varB3"))
  )

  expected <- matrix(
    data = c(
      perfMix["mix2"],
      perfPur["varC2"] * mix2pur[["mix2"]]["varC2"] +
        perfPur["varA1"] * mix2pur[["mix2"]]["varA1"] +
        perfPur["varB3"] * mix2pur[["mix2"]]["varB3"],
      perfMix["mix2"] /
        (perfPur["varC2"] * mix2pur[["mix2"]]["varC2"] +
          perfPur["varA1"] * mix2pur[["mix2"]]["varA1"] +
          perfPur["varB3"] * mix2pur[["mix2"]]["varB3"]),
      perfMix["mix1"],
      perfPur["varB3"] * mix2pur[["mix1"]]["varB3"] +
        perfPur["varC2"] * mix2pur[["mix1"]]["varC2"],
      perfMix["mix1"] /
        (perfPur["varB3"] * mix2pur[["mix1"]]["varB3"] +
          perfPur["varC2"] * mix2pur[["mix1"]]["varC2"])
    ),
    nrow = 2, ncol = 3, byrow = TRUE,
    dimnames = list(
      c("mix2", "mix1"),
      c("mix", "weightedPurs", "RY_RD18")
    )
  )

  observed <- RY_RD18(mixYields = perfMix, monoYields = perfPur, mix2pur)

  expect_equal(observed, expected)
})

test_that("overyielding", {
  perfMix <- setNames(
    c(40, 50),
    c("mix2", "mix1")
  )
  perfPur <- setNames(
    c(70, 35, 20),
    c("varA1", "varC2", "varB3")
  )
  mix2pur <- list(
    "mix1" = c("varB3", "varC2"),
    "mix2" = c("varC2", "varA1", "varB3")
  )

  expected <- matrix(
    data = c(
      perfMix["mix2"],
      mean(perfPur[mix2pur[["mix2"]]]),
      perfMix["mix2"] / mean(perfPur[mix2pur[["mix2"]]]),
      perfMix["mix1"],
      mean(perfPur[mix2pur[["mix1"]]]),
      perfMix["mix1"] / mean(perfPur[mix2pur[["mix1"]]])
    ),
    nrow = 2, ncol = 3, byrow = TRUE,
    dimnames = list(
      c("mix2", "mix1"),
      c("mix", "meanPur", "OY")
    )
  )

  observed <- plantmix:::overyielding(perfMix = perfMix, perfPur = perfPur, mix2pur)

  expect_equal(observed, expected)
})

test_that("LER", {
  dat <- data.frame(
    solecrop = c(5, 15),
    intercrop = c(4, 9),
    row.names = c("grain", "fruit")
  )

  expected <- list(
    pLER = setNames(c(0.8, 0.6), c("grain", "fruit")),
    LER = 1.4
  )
  observed <- LER(dat)

  expect_equal(observed, expected)
})

test_that("CC", {
  ## example of a binary varietal mixture with unequal proportions:
  dat <- data.frame(
    ID = c("geno1", "geno2", "mixg1g2", "mixg1g2"),
    focal = c("geno1", "geno2", "geno1", "geno2"),
    prop = c(1, 1, 0.6, 0.4),
    yield = c(50, 40, 25, 22)
  )
  expected <- data.frame(dat,
    CC = c(
      NA,
      NA,
      (25 / 22) / ((0.6 * 50) / (0.4 * 40)) - 1,
      (22 / 25) / ((0.4 * 40) / (0.6 * 50)) - 1
    )
  )
  observed <- CC(dat)
  expect_equal(observed, expected)
})

test_that("RII", {
  ## example of a binary varietal mixture
  sow_density <- 200
  dat <- data.frame(
    ID = c("geno1", "geno2", "mixg1g2", "mixg1g2"),
    focal = c("geno1", "geno2", "geno1", "geno2"),
    prop = c(1, 1, 0.6, 0.4),
    yield_qt_ha = c(50, 40, 25, 22)
  )
  dat$yield_g_m2 <- (dat$yield_qt_ha / 10^4) * 10^6
  dat$yield_g_plant <- dat$yield_g_m2 / (sow_density * dat$prop)

  expected <- data.frame(dat,
    RII = c(1, 1, NA, NA)
  )
  i <- 0
  Ymix <- dat$yield_g_plant[3 + i]
  Ymono <- dat$yield_g_plant[1 + i]
  expected$RII[3 + i] <- (Ymix - Ymono) / (Ymix + Ymono)
  i <- 1
  Ymix <- dat$yield_g_plant[3 + i]
  Ymono <- dat$yield_g_plant[1 + i]
  expected$RII[3 + i] <- (Ymix - Ymono) / (Ymix + Ymono)

  observed <- RII(dat, colY = "yield_g_plant")

  expect_equal(observed, expected)
})

test_that("RIInet", {
  ## example of a binary varietal mixture
  sow_density <- 200
  dat <- data.frame(
    ID = c("geno1", "geno2", "mixg1g2", "mixg1g2"),
    focal = c("geno1", "geno2", "geno1", "geno2"),
    prop = c(1, 1, 0.6, 0.4),
    yield_qt_ha = c(50, 40, 25, 22)
  )
  dat$yield_g_m2 <- (dat$yield_qt_ha / 10^4) * 10^6
  dat$yield_g_plant <- dat$yield_g_m2 / (sow_density * dat$prop)

  dat <- RII(dat, colY = "yield_g_plant")

  expected <- data.frame(
    ID = "mixg1g2",
    RIInet = dat$prop[3] * dat$RII[3] + dat$prop[4] * dat$RII[4]
  )

  observed <- RIInet(dat)

  expect_equal(observed, expected)
})

test_that("ratioComp", {
  # without names
  yieldPure <- c(12, 16)
  yieldMix <- c(14, 15)
  expected <- c(14 * 2 / 12, 15 * 2 / 16)
  observed <- ratioComp(yieldPure, yieldMix)
  expect_equal(expected, observed)

  # with names
  yieldPure <- c("var37" = 12, "var12" = 16)
  yieldMix <- c("var12" = 15, "var37" = 14)
  expected <- c("var37" = 14 * 2 / 12, "var12" = 15 * 2 / 16)
  observed <- ratioComp(yieldPure, yieldMix)
  expect_equal(expected, observed)
})

test_that("ratioProd", {
  input.vec1 <- c(12, 16)
  expected <- (12 + 16) / 2
  observed <- ratioProd(input.vec1)
  expect_equal(expected, observed)
})

test_that("ratioAggr", {
  input.vec1 <- c(12, 16)
  expected <- c(12 - 16, 16 - 12)
  observed <- ratioAggr(input.vec1)
  expect_equal(expected, observed)
})

test_that("ratioRelAggr", {
  input.vec1 <- c(12, 16)
  expected <- c(12 / 16, 16 / 12)
  observed <- ratioRelAggr(input.vec1)
  expect_equal(expected, observed)
})
