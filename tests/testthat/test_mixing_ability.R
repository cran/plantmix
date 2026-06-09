library(plantmix)
context("mixing_ability")

## ===========================================================================
## Tests of functions useful to preprocess the input data

test_that("getMixedGenos_sep-comma", {
  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3,var2", "var1,var3", "var1,var2"),
    yield = c(10, 7, 9)
  )
  expected <- c("var1", "var2", "var3")
  observed <- getMixedGenos(df = dat, col = "varieties", sep = ",")
  expect_equal(observed, expected)
})

test_that("getMixedGenos_sep-dash", {
  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3-var2", "var1-var3", "var1-var2"),
    yield = c(10, 7, 9)
  )
  expected <- c("var1", "var2", "var3")
  observed <- getMixedGenos(df = dat, col = "varieties", sep = "-")
  expect_equal(observed, expected)
})

test_that("getGenosPerMix_sep-comma", {
  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3,var2", "var1,var3", "var1,var2"),
    yield = c(10, 7, 9)
  )
  expected <- list(
    c("var2", "var3"),
    c("var1", "var3"),
    c("var1", "var2")
  )
  observed <- getGenosPerMix(df = dat, col = "varieties", sep = ",")
  expect_equal(observed, expected)
})

test_that("getGenosPerMix_sep-underscore", {
  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3_var2", "var1_var3", "var1_var2"),
    yield = c(10, 7, 9)
  )
  expected <- list(
    c("var2", "var3"),
    c("var1", "var3"),
    c("var1", "var2")
  )
  observed <- getGenosPerMix(df = dat, col = "varieties", sep = "_")
  expect_equal(observed, expected)
})

test_that("getGenoPairs_binary-mixtures", {
  mixed.genos <- c("var1", "var2", "var3", "var10")
  genos.per.mix <- list(
    c("var2", "var3"),
    c("var1", "var3"),
    c("var2", "var10")
  )
  expected <- c(
    "var1-var1", "var1-var3",
    "var10-var10", "var2-var10",
    "var2-var2", "var2-var3",
    "var3-var3"
  )
  observed <- getGenoPairs(
    mixed.genos = mixed.genos,
    genos.per.mix = genos.per.mix,
    sep = "-"
  )
  expect_equal(observed, expected)
})

test_that("getGenoPairs_binary-ternary-mixtures", {
  mixed.genos <- c("var1", "var2", "var3", "var10")
  genos.per.mix <- list(
    c("var2", "var3"),
    c("var1", "var3", "var10")
  )
  expected <- c(
    "var1-var1", "var1-var10", "var1-var3",
    "var10-var10",
    "var2-var2", "var2-var3",
    "var3-var10", "var3-var3"
  )
  observed <- getGenoPairs(
    mixed.genos = mixed.genos,
    genos.per.mix = genos.per.mix,
    sep = "-"
  )
  expect_equal(observed, expected)
})

test_that("getPhenos_one-trait", {
  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3_var2", "var1_var3", "var1_var2"),
    yield = c(10, 7, 9)
  )
  expected <- list("yield" = c(10, 7, 9))
  observed <- getPhenos(df = dat, response.names = "yield")
  expect_equal(observed, expected)
})

test_that("getPhenos_two-traits", {
  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3_var2", "var1_var3", "var1_var2"),
    yield = c(10, 7, 9),
    disease = c(0.1, 0.5, 0.2)
  )
  expected <- list(
    "yield" = c(10, 7, 9),
    "disease" = c(0.1, 0.5, 0.2)
  )
  observed <- getPhenos(df = dat, response.names = c("yield", "disease"))
  expect_equal(observed, expected)
})

test_that("getMixtureList_vec", {
  mixtures <- c(
    "comp1+comp2",
    "comp1+comp3"
  )
  expected <- list(
    "comp1+comp2" = c("comp1", "comp2"),
    "comp1+comp3" = c("comp1", "comp3")
  )
  observed <- getMixtureList(mixtures, "+")
  expect_equal(observed, expected)
})

test_that("getMixtureList_mat", {
  mixtures <- matrix(
    c(
      "comp1", "comp2",
      "comp1", "comp3"
    ),
    nrow = 2, ncol = 2, byrow = TRUE,
    dimnames = list(
      c("mix1", "mix2"),
      c("comp1", "comp2")
    )
  )
  expected <- list(
    "mix1" = c("comp1", "comp2"),
    "mix2" = c("comp1", "comp3")
  )
  observed <- getMixtureList(mixtures)
  expect_equal(observed, expected)
})

test_that("getMixtureList_df", {
  mixtures <- data.frame(
    "comp1" = c("comp1", "comp3"),
    "comp2" = c("comp2", "comp1")
  )
  rownames(mixtures) <- c("mix1", "mix2")
  expected <- list(
    "mix1" = c("comp1", "comp2"),
    "mix2" = c("comp3", "comp1")
  )
  observed <- getMixtureList(mixtures)
  expect_equal(observed, expected)
})

test_that("getMixtureList_list", {
  mixtures <- list(
    "mix1" = c("comp1", "comp2"),
    "mix2" = c("comp1", "comp3")
  )
  expected <- mixtures
  observed <- getMixtureList(mixtures)
  expect_equal(observed, expected)
})

test_that("getMixtureList_list_fact", {
  mixtures <- list(
    "mix1" = factor(c("comp1", "comp2")),
    "mix2" = factor(c("comp1", "comp3"))
  )
  expected <- list(
    "mix1" = c("comp1", "comp2"),
    "mix2" = c("comp1", "comp3")
  )
  observed <- getMixtureList(mixtures)
  expect_equal(observed, expected)
})

test_that("getMixturesPerGeno", {
  mix2genos <- list(
    "mix1" = c("geno1", "geno2"),
    "mix2" = c("geno1", "geno3")
  )
  expected <- list(
    "geno1" = c("mix1", "mix2"),
    "geno2" = "mix1",
    "geno3" = "mix2"
  )
  observed <- getMixturesPerGeno(mix2genos)
  expect_equal(observed, expected)
})

test_that("pivotMixData2Long_speciesMix_pureAndBinMix_noYield", {
  dat <- data.frame(
    name = c("wheat01-pea02", "wheat01-pea01", "wheat01"),
    stand = c("inter", "inter", "sole"),
    x = c(1, 1, 1),
    y = c(1, 2, 3),
    check.names = FALSE,
    stringsAsFactors = TRUE
  )

  expected <- data.frame(
    focal = c("wheat01", "pea02", "wheat01", "pea01", "wheat01"),
    neighbor = c("pea02", "wheat01", "pea01", "wheat01", "wheat01"),
    name = c(
      rep("wheat01-pea02", 2),
      rep("wheat01-pea01", 2),
      "wheat01"
    ),
    stand = c(rep("inter", 4), "sole"),
    x = 1,
    y = c(1, 1, 2, 2, 3),
    stringsAsFactors = TRUE
  )

  observed <- pivotMixData2Long(
    df = dat, colC = "name", sep = "-",
    genos = list(
      "cereal" = c("wheat01"),
      "legume" = c("pea01", "pea02")
    )
  )

  expect_equal(observed, expected)
})

test_that("pivotMixData2Long_speciesMix_pureAndBinMix", {
  dat <- data.frame(
    name = c("wheat01-pea02", "wheat01-pea01", "wheat01"),
    stand = c("inter", "inter", "sole"),
    x = c(1, 1, 1),
    y = c(1, 2, 3),
    "yield_cereal" = c(10, 11, 20),
    "yield_legume" = c(9, 8, NA),
    check.names = FALSE,
    stringsAsFactors = TRUE
  )

  expected <- data.frame(
    focal = c("wheat01", "pea02", "wheat01", "pea01", "wheat01"),
    neighbor = c("pea02", "wheat01", "pea01", "wheat01", "wheat01"),
    name = c(
      rep("wheat01-pea02", 2),
      rep("wheat01-pea01", 2),
      "wheat01"
    ),
    stand = c(rep("inter", 4), "sole"),
    x = 1,
    y = c(1, 1, 2, 2, 3),
    yield = c(10, 9, 11, 8, 20),
    stringsAsFactors = TRUE
  )

  observed <- pivotMixData2Long(
    df = dat, colC = "name", sep = "-", prefixY = "yield", sepY = "_",
    genos = list(
      "cereal" = c("wheat01"),
      "legume" = c("pea01", "pea02")
    )
  )

  expect_equal(observed, expected)
})

test_that("pivotMixData2Long_varMix_pureAndBinMix", {
  dat <- data.frame(
    name = c("g1-g2", "g1-g2", "g1", "g2"),
    stand = c("mix", "mix", "mono", "mono"),
    block = c("A", "B", "A", "A"),
    x = c(1, 2, 1, 1),
    y = c(1, 1, 2, 3),
    "yield_focal" = c(10, 12, 20, 18),
    "yield_neighbor" = c(11, 9, NA, NA),
    check.names = FALSE,
    stringsAsFactors = TRUE
  )

  expected <- data.frame(
    focal = c("g1", "g2", "g1", "g2", "g1", "g2"),
    neighbor = c("g2", "g1", "g2", "g1", "g1", "g2"),
    name = c(
      rep("g1-g2", 2),
      rep("g1-g2", 2),
      "g1", "g2"
    ),
    stand = c(rep("mix", 4), rep("mono", 2)),
    block = c(rep("A", 2), rep("B", 2), rep("A", 2)),
    x = c(1, 1, 2, 2, 1, 1),
    y = c(1, 1, 1, 1, 2, 3),
    yield = c(10, 11, 12, 9, 20, 18),
    stringsAsFactors = TRUE
  )

  observed <- pivotMixData2Long(
    df = dat, colC = "name", sep = "-", prefixY = "yield", sepY = "_",
    genos = list(c("g1", "g2"))
  )

  expect_equal(observed, expected)
})

test_that("pivotMixData2Wide_binMix-only", {
  dat <- data.frame(
    focal = c("wheat01", "pea02", "wheat01", "pea01"),
    neighbor = c("pea02", "wheat01", "pea01", "wheat01"),
    name = c(
      rep("wheat01-pea02", 2),
      rep("wheat01-pea01", 2)
    ),
    stand = rep("inter", 4),
    block = rep("A", 4),
    x = 1,
    y = c(1, 1, 2, 2),
    yield = c(10, 11, 12, 13),
    stringsAsFactors = TRUE
  )

  expected <- data.frame(
    focal = c("wheat01", "wheat01"),
    neighbor = c("pea02", "pea01"),
    name = c("wheat01-pea02", "wheat01-pea01"),
    stand = c("inter", "inter"),
    block = c("A", "A"),
    x = c(1, 1),
    y = c(1, 2),
    "yield_focal" = c(10, 12),
    "yield_neighbor" = c(11, 13),
    check.names = FALSE,
    stringsAsFactors = TRUE
  )

  observed <- pivotMixData2Wide(
    df = dat, colIDstand = "name",
    colIDfocal = "focal",
    colIDneighbors = "neighbor",
    colPlot = c("x", "y"),
    colY = "yield", sepY = "_"
  )

  expect_equal(observed, expected)
})

test_that("pivotMixData2Wide_binMix-only_no-neighbor-col", {
  dat <- data.frame(
    geno = c("wheat01", "pea02", "wheat01", "pea01"),
    standID = c(
      rep("wheat01+pea02", 2),
      rep("wheat01+pea01", 2)
    ),
    stand = rep("inter", 4),
    block = rep("A", 4),
    x = 1,
    y = c(1, 1, 2, 2),
    yield = c(10, 11, 12, 13),
    stringsAsFactors = TRUE
  )

  expected <- data.frame(
    geno = c("wheat01", "wheat01"),
    standID = c("wheat01+pea02", "wheat01+pea01"),
    stand = c("inter", "inter"),
    block = c("A", "A"),
    x = c(1, 1),
    y = c(1, 2),
    "neighbor" = c("pea02", "pea01"),
    "yield_geno" = c(10, 12),
    "yield_neighbor" = c(11, 13),
    check.names = FALSE,
    stringsAsFactors = TRUE
  )

  observed <- pivotMixData2Wide(
    df = dat, colIDstand = "standID",
    colIDfocal = "geno",
    colIDneighbors = "neighbor", sepFocalNeighbors = "+",
    colPlot = c("x", "y"),
    colY = "yield", sepY = "_"
  )

  expect_equal(observed, expected)
})

test_that("pivotMixData2Wide_monoAndBinMix", {
  dat <- data.frame(
    focal = c("wheat01", "pea02", "wheat01", "pea01", "wheat01"),
    neighbor = c("pea02", "wheat01", "pea01", "wheat01", "wheat01"),
    name = c(
      rep("wheat01-pea02", 2),
      rep("wheat01-pea01", 2),
      "wheat01"
    ),
    stand = c(rep("inter", 4), "sole"),
    block = "A",
    x = 1,
    y = c(1, 1, 2, 2, 3),
    yield = c(10, 11, 12, 13, 20.5),
    stringsAsFactors = TRUE
  )

  expected <- data.frame(
    focal = c("wheat01", "wheat01", "wheat01"),
    neighbor = c("pea02", "pea01", "wheat01"),
    name = c("wheat01-pea02", "wheat01-pea01", "wheat01"),
    stand = c("inter", "inter", "sole"),
    block = "A",
    x = 1,
    y = c(1, 2, 3),
    yield_focal = c(10, 12, 20.5),
    yield_neighbor = c(11, 13, NA),
    stringsAsFactors = TRUE
  )

  observed <- pivotMixData2Wide(
    df = dat, colIDstand = "name",
    colIDfocal = "focal",
    colIDneighbors = "neighbor",
    colPlot = c("x", "y"),
    colY = "yield", sepY = "_"
  )

  expect_equal(observed, expected)
})

## ===========================================================================
## Tests of the function making the design matrix for the GMAs

test_that("mkZGMA_only2wayMixtures", {
  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3-var2", "var1-var3", "var1-var2"),
    pheno = c(10, 7, 9)
  )
  expected <- matrix(
    c(
      0, 0.5, 0.5,
      0.5, 0, 0.5,
      0.5, 0.5, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 3
  )
  colnames(expected) <- paste0("var", 1:3)

  ## Tim's code
  observed <- mkZGMA(df = dat, col = "varieties", sep = "-")
  expect_equal(observed, expected)

  ## Emma's code
  genos <- getMixedGenos(df = dat, col = "varieties", sep = "-")
  observed <- BuildZ_GMA(
    df = dat, Gcol = which(colnames(dat) == "varieties"),
    G = genos, sep = "-"
  )
  expect_equal(observed, expected)
})

test_that("mkZGMA_2and3wayMixtures", {
  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3-var2", "var1-var3", "var1-var2-var3"),
    pheno = c(10, 7, 9)
  )
  expected <- matrix(
    c(
      0, 0.5, 0.5,
      0.5, 0, 0.5,
      1 / 3, 1 / 3, 1 / 3
    ),
    byrow = TRUE, nrow = 3, ncol = 3
  )
  colnames(expected) <- paste0("var", 1:3)

  ## Tim's code
  observed <- mkZGMA(df = dat, col = "varieties", sep = "-")
  expect_equal(observed, expected)

  ## Emma's code
  genos <- getMixedGenos(df = dat, col = "varieties", sep = "-")
  observed <- BuildZ_GMA(
    df = dat, Gcol = which(colnames(dat) == "varieties"),
    G = genos, sep = "-"
  )
  expect_equal(observed, expected)
})

test_that("mkZGMA_pureStandsAnd2wayMixtures", {
  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3", "var1-var3", "var1-var2"),
    pheno = c(10, 7, 9)
  )
  expected <- matrix(
    c(
      0, 0, 1,
      0.5, 0, 0.5,
      0.5, 0.5, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 3
  )
  colnames(expected) <- paste0("var", 1:3)

  ## Tim's code
  observed <- mkZGMA(df = dat, col = "varieties", sep = "-")
  expect_equal(observed, expected)

  ## Emma's code
  genos <- getMixedGenos(df = dat, col = "varieties", sep = "-")
  observed <- BuildZ_GMA(
    df = dat, Gcol = which(colnames(dat) == "varieties"),
    G = genos, sep = "-"
  )
  expect_equal(observed, expected)
})

## ===========================================================================
## Tests of the function making the design matrix for the SMAs

## -------------------------------------------------------------------
## input = only 2-way mixtures

test_that("mkZSMA_only2wayMixtures_with-SMAij_without-SMAii", {
  ## case corresponding to the classical model

  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3-var2", "var1-var3", "var1-var2"),
    pheno = c(10, 7, 9)
  )

  expected <- matrix(
    c(
      0, 0, 1,
      0, 1, 0,
      1, 0, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 3
  )
  colnames(expected) <- c("var1-var2", "var1-var3", "var2-var3")

  ## Tim's code
  observed <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "no"
  )
  expect_equal(observed, expected)

  ## Emma's code: it cannot do that as it always includes columns
  ## for pure stands
})

test_that("mkZSMA_only2wayMixtures_with-SMAij_with-SMAii-onlyPur", {
  ## case corresponding to model 2 from Forst et al (2019)

  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3-var2", "var1-var3", "var1-var2-var3"),
    pheno = c(10, 7, 9)
  )

  ## Tim's code: note that pure stands are NOT included as it is model 2
  expected <- matrix(
    c(
      0, 0, 1,
      0, 1, 0,
      1 / 3, 1 / 3, 1 / 3
    ),
    byrow = TRUE, nrow = 3, ncol = 3
  )
  colnames(expected) <- c("var1-var2", "var1-var3", "var2-var3")
  observed <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "only_pur"
  )
  expect_equal(observed, expected)

  ## Emma's code: note that pure stands are included even if it is model 2
  expected <- matrix(
    c(
      0, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0,
      0, 1 / 3, 1 / 3, 0, 1 / 3, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 6
  )
  colnames(expected) <- c(
    "var1-var1", "var1-var2", "var1-var3",
    "var2-var2", "var2-var3", "var3-var3"
  )
  genos <- getMixedGenos(df = dat, col = "varieties", sep = "-")
  genos_per_mix <- getGenosPerMix(dat, "varieties", "-")
  geno_pairs <- getGenoPairs(genos, genos_per_mix)
  observed <- BuildZ_SMAmodel2(
    df = dat, Gcol = which(colnames(dat) == "varieties"),
    G = genos, Gpairs = geno_pairs, sep = "-"
  )
  expect_equal(observed, expected)
})

test_that("mkZSMA_only2wayMixtures_with-SMAij_with-SMAii-pureAndMix", {
  ## case corresponding to model 3 from Forst et al (2019)

  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3-var2", "var1-var3", "var1-var2"),
    pheno = c(10, 7, 9)
  )

  expected <- matrix(
    c(
      0, 0, 0, 0.25, 0.5, 0.25,
      0.25, 0, 0.5, 0, 0, 0.25,
      0.25, 0.5, 0, 0.25, 0, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 6
  )
  colnames(expected) <- c(
    "var1-var1", "var1-var2", "var1-var3",
    "var2-var2", "var2-var3", "var3-var3"
  )

  ## Tim's code
  observed <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "pur_mix"
  )
  expect_equal(observed, expected)

  ## Emma's code
  genos <- plantmix:::getMixedGenos(df = dat, col = "varieties", sep = "-")
  genos_per_mix <- plantmix:::getGenosPerMix(dat, "varieties", "-")
  geno_pairs <- plantmix:::getGenoPairs(genos, genos_per_mix)
  observed <- plantmix:::BuildZ_SMAmodel3(
    df = dat, Gcol = which(colnames(dat) == "varieties"),
    Gpairs = geno_pairs, sep = "-"
  )
  expect_equal(observed, expected)
})

## -------------------------------------------------------------------
## input = 2-way and 3-way mixtures

test_that("mkZSMA_2and3wayMixtures_with-SMAij_without-SMAii", {
  ## case corresponding to the classical model

  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3-var2", "var1-var3", "var1-var2-var3"),
    pheno = c(10, 7, 9)
  )

  expected <- matrix(
    c(
      0, 0, 1,
      0, 1, 0,
      1 / 3, 1 / 3, 1 / 3
    ),
    byrow = TRUE, nrow = 3, ncol = 3
  )
  colnames(expected) <- c("var1-var2", "var1-var3", "var2-var3")

  ## Tim's code
  observed <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "no"
  )
  expect_equal(observed, expected)

  ## Emma's code: it cannot do that as it always includes columns
  ## for pure stands
})

test_that("mkZSMA_2and3wayMixtures_with-SMAij_with-SMAii-onlyPur", {
  ## case corresponding to model 2 from Forst et al (2019)

  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3-var2", "var1-var3", "var1-var2-var3"),
    pheno = c(10, 7, 9)
  )

  expected <- matrix(
    c(
      0, 0, 1,
      0, 1, 0,
      1 / 3, 1 / 3, 1 / 3
    ),
    byrow = TRUE, nrow = 3, ncol = 3
  )
  colnames(expected) <- c("var1-var2", "var1-var3", "var2-var3")
  ## this specific data set has no pure stand, hence the design matrix Z_SMA
  ## has no column for them even though this model can include SMA_ii in them

  ## Tim's code
  observed <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "only_pur"
  )
  expect_equal(observed, expected)

  ## Emma's code
  genos <- plantmix:::getMixedGenos(df = dat, col = "varieties", sep = "-")
  genos_per_mix <- plantmix:::getGenosPerMix(dat, "varieties", "-")
  geno_pairs <- plantmix:::getGenoPairs(genos, genos_per_mix)
  observed <- plantmix:::BuildZ_SMAmodel2(
    df = dat, Gcol = which(colnames(dat) == "varieties"),
    G = genos, Gpairs = geno_pairs, sep = "-"
  )
  isEmpty <- (colSums(observed) == 0)
  observed <- observed[, !isEmpty]
  expect_equal(observed, expected)
})

test_that("mkZSMA_2and3wayMixtures_with-SMAij_with-SMAii-pureAndMix", {
  ## case corresponding to model 3 from Forst et al (2019)

  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3-var2", "var1-var3", "var1-var2-var3"),
    pheno = c(10, 7, 9)
  )

  expected <- matrix(
    c(
      0, 0, 0, 0.25, 0.5, 0.25,
      0.25, 0, 0.5, 0, 0, 0.25,
      1 / 3^2, 2 / 3^2, 2 / 3^2, 1 / 3^2, 2 / 3^2, 1 / 3^2
    ),
    byrow = TRUE, nrow = 3, ncol = 6
  )
  colnames(expected) <- c(
    "var1-var1", "var1-var2", "var1-var3",
    "var2-var2", "var2-var3", "var3-var3"
  )

  ## Tim's code
  observed <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "pur_mix"
  )
  expect_equal(observed, expected)

  ## Emma's code
  genos <- plantmix:::getMixedGenos(df = dat, col = "varieties", sep = "-")
  genos_per_mix <- plantmix:::getGenosPerMix(dat, "varieties", "-")
  geno_pairs <- plantmix:::getGenoPairs(genos, genos_per_mix)
  observed <- plantmix:::BuildZ_SMAmodel3(
    df = dat, Gcol = which(colnames(dat) == "varieties"),
    Gpairs = geno_pairs, sep = "-"
  )
  expect_equal(observed, expected)
})

## -------------------------------------------------------------------
## input = 2-way and 3-way mixtures as well as pure stands

test_that("mkZSMA_pureStandsAnd2wayMixtures_with-SMAij_without-SMAii", {
  ## case corresponding to the classical model

  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3", "var1-var3", "var1-var2"),
    pheno = c(10, 7, 9)
  )

  expected <- matrix(
    c(
      0, 0,
      0, 1,
      1, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 2
  )
  colnames(expected) <- c("var1-var2", "var1-var3")

  ## Tim's code: skipping unused columns
  observed <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "no"
  )
  expect_equal(observed, expected)

  ## Tim's code: keeping unused columns
  expected2 <- matrix(
    c(
      0, 0, 0,
      0, 1, 0,
      1, 0, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 3
  )
  colnames(expected2) <- c("var1-var2", "var1-var3", "var2-var3")
  observed2 <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "no",
    skipUnusedCols = FALSE
  )
  expect_equal(observed2, expected2)

  ## Emma's code: it cannot do that as it always includes columns
  ## for pure stands
})

test_that("mkZSMA_pureStandsAnd2wayMixtures_with-SMAij_with-SMAii-onlyPur", {
  ## case corresponding to model 2 from Forst et al (2019)

  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3", "var1-var3", "var1-var2"),
    pheno = c(10, 7, 9)
  )

  expected <- matrix(
    c(
      0, 0, 1,
      0, 1, 0,
      1, 0, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 3
  )
  colnames(expected) <- c("var1-var2", "var1-var3", "var3-var3")

  ## Tim's code
  observed <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "only_pur"
  )
  expect_equal(observed, expected)

  ## Emma's code
  genos <- plantmix:::getMixedGenos(df = dat, col = "varieties", sep = "-")
  genos_per_mix <- plantmix:::getGenosPerMix(dat, "varieties", "-")
  geno_pairs <- plantmix:::getGenoPairs(genos, genos_per_mix)
  observed <- plantmix:::BuildZ_SMAmodel2(
    df = dat, Gcol = which(colnames(dat) == "varieties"),
    G = genos, Gpairs = geno_pairs, sep = "-"
  )
  isEmpty <- (colSums(observed) == 0)
  observed <- observed[, !isEmpty]
  expect_equal(observed, expected)
})

test_that("mkZSMA_pureStandsAnd2wayMixtures_with-SMAij_with-SMAii-pureAndMix", {
  ## case corresponding to model 3 from Forst et al (2019)

  dat <- data.frame(
    mix = paste0("mix", 1:3),
    varieties = c("var3", "var1-var3", "var1-var2"),
    pheno = c(10, 7, 9)
  )

  expected <- matrix(
    c(
      0, 0, 0, 0, 1,
      0.25, 0, 0.5, 0, 0.25,
      0.25, 0.5, 0, 0.25, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 5
  )
  colnames(expected) <- c(
    "var1-var1", "var1-var2", "var1-var3",
    "var2-var2", "var3-var3"
  )

  ## Tim's code: skipping unused columns
  observed <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "pur_mix"
  )
  expect_equal(observed, expected)

  ## Emma's code
  genos <- plantmix:::getMixedGenos(df = dat, col = "varieties", sep = "-")
  genos_per_mix <- plantmix:::getGenosPerMix(dat, "varieties", "-")
  geno_pairs <- plantmix:::getGenoPairs(genos, genos_per_mix)
  observed <- plantmix:::BuildZ_SMAmodel3(
    df = dat, Gcol = which(colnames(dat) == "varieties"),
    Gpairs = geno_pairs, sep = "-"
  )
  expect_equal(observed, expected)

  ## Tim's code: keeping unused columns
  expected2 <- matrix(
    c(
      0, 0, 0, 0, 0, 1,
      0.25, 0, 0.5, 0, 0, 0.25,
      0.25, 0.5, 0, 0.25, 0, 0
    ),
    byrow = TRUE, nrow = 3, ncol = 6
  )
  colnames(expected2) <- c(
    "var1-var1", "var1-var2", "var1-var3",
    "var2-var2", "var2-var3", "var3-var3"
  )
  observed2 <- mkZSMA(
    df = dat, col = "varieties", sep = "-",
    inc_SMA_ii = "pur_mix",
    skipUnusedCols = FALSE
  )
  expect_equal(observed2, expected2)
})

## ===========================================================================
## Tests of the function making the design matrices for DBVs and SBVs

test_that("mkZinterSpe", {
  dat <- data.frame(
    focal = c(
      "wheat01", "pea02",
      "wheat01", "pea01",
      "wheat01"
    ),
    neighbors = c(
      "pea02", "wheat01",
      "pea01", "wheat01",
      "wheat01"
    ),
    name = c(
      "wheat01-pea02", "wheat01-pea02",
      "wheat01-pea01", "wheat01-pea01",
      "wheat01-wheat01"
    ),
    stand = c("inter", "inter", "inter", "inter", "sole"),
    species = c("wheat", "pea", "wheat", "pea", "wheat")
  )

  ## expected wheat matrices:
  exp_wheat_Z_F <- matrix(
    c(
      1,
      0,
      1,
      0,
      1
    ),
    nrow = 5, ncol = 1, byrow = TRUE,
    dimnames = list(
      as.character(1:5),
      "wheat01"
    )
  )
  exp_wheat_Z_N <- matrix(
    c(
      0,
      1,
      0,
      1,
      0
    ), # <- because it is a sole crop, see Z_SIS
    nrow = 5, ncol = 1, byrow = TRUE,
    dimnames = list(
      as.character(1:5),
      "wheat01"
    )
  )
  exp_wheat <- list("Z_DBV" = exp_wheat_Z_F, "Z_SBV" = exp_wheat_Z_N)

  ## expected pea matrices:
  exp_pea_Z_F <- matrix(
    c(
      0, 0,
      0, 1,
      0, 0,
      1, 0,
      0, 0
    ),
    nrow = 5, ncol = 2, byrow = TRUE,
    dimnames = list(
      as.character(1:5),
      c("pea01", "pea02")
    )
  )
  exp_pea_Z_N <- matrix(
    c(
      0, 1,
      0, 0,
      1, 0,
      0, 0,
      0, 0
    ),
    nrow = 5, ncol = 2, byrow = TRUE,
    dimnames = list(
      as.character(1:5),
      c("pea01", "pea02")
    )
  )
  exp_pea <- list("Z_DBV" = exp_pea_Z_F, "Z_SBV" = exp_pea_Z_N)

  expected <- list("wheat" = exp_wheat, "pea" = exp_pea)

  observed <- mkZinterspe(dat,
    genosPerSp = list(
      "wheat" = "wheat01",
      "pea" = c("pea01", "pea02")
    ),
    colIDfocal = "focal", colIDneighbors = "neighbors"
  )

  expect_equal(observed, expected)
})

test_that("mkZSIS", {
  dat <- data.frame(
    focal = c(
      "wheat01", "pea02",
      "wheat01", "pea01",
      "wheat01"
    ),
    neighbors = c(
      "pea02", "wheat01",
      "pea01", "wheat01",
      "wheat01"
    ),
    name = c(
      "wheat01-pea02", "wheat01-pea02",
      "wheat01-pea01", "wheat01-pea01",
      "wheat01-wheat01"
    ),
    stand = c("inter", "inter", "inter", "inter", "sole"),
    species = c("wheat", "pea", "wheat", "pea", "wheat")
  )

  ## Z_SIS_wheat
  expected <- matrix(
    c(
      0.5,
      0,
      0.5,
      0,
      1
    ),
    nrow = 5, ncol = 1, byrow = TRUE
  )
  colnames(expected) <- "wheat01"
  observed <- mkZSIS(dat,
    colIDfocal = "focal", colStand = "stand",
    colS = "species", species = "wheat",
    weight = 0.5
  )
  expect_equal(observed, expected)

  ## Z_SIS_pea
  expected <- matrix(
    c(
      0, 0,
      0, 0.25,
      0, 0,
      0.25, 0,
      0, 0
    ),
    nrow = 5, ncol = 2, byrow = TRUE
  )
  colnames(expected) <- c("pea01", "pea02")
  observed <- mkZSIS(dat,
    colIDfocal = "focal", colStand = "stand",
    colS = "species", species = "pea",
    weight = 0.25
  )
  expect_equal(observed, expected)
})

## ===========================================================================
## Tests of the wrapping of lme4 functions for GMA-SMA models

simulDat <- function() {
  ## simulate data: y = X beta + Z GMA + e
  nbGenos <- 25
  genos <- sprintf("g%02i", 1:nbGenos)
  pairs <- t(combn(x = genos, m = 2))
  stands <- paste(pairs[, 1], pairs[, 2], sep = "_")
  nbBlocks <- 3
  blocks <- LETTERS[1:nbBlocks]
  dat <- do.call(rbind, lapply(blocks, function(block) {
    cbind(
      stands = as.data.frame(stands, stringsAsFactors = TRUE),
      block = as.factor(block)
    )
  }))
  listContr <- list(block = "contr.sum")
  X <- model.matrix(~ 1 + block, data = dat, contrasts = listContr)
  Z_GMA <- mkZGMA(df = dat, col = "stands", sep = "_")
  truth <- list(
    "intercept" = 100,
    "var_GMA" = 10,
    "var_error" = 1
  )
  set.seed(1234)
  truth[["blockEffs"]] <- sample(x = c(-1, 1), size = nbBlocks - 1, replace = TRUE) *
    rnorm(n = nbBlocks - 1, mean = 3, sd = 5)
  truth[["GMAs"]] <- rnorm(n = nbGenos, mean = 0, sd = sqrt(truth$var_GMA))
  truth[["errors"]] <- rnorm(n = nrow(dat), mean = 0, sd = sqrt(truth$var_error))
  y <- X %*% c(truth$intercept, truth$blockEffs) +
    Z_GMA %*% truth$GMAs +
    truth$errors
  dat$pheno <- y[, 1]
  if (FALSE) {
    str(dat)
    hist(dat$pheno, las = 1)
    boxplot(pheno ~ block, data = dat, las = 1)
  }
  return(list(
    "truth" = truth, "dat" = dat, "listContr" = listContr,
    "X" = X, "Z_GMA" = Z_GMA
  ))
}

## test_that("lmerGMASMA_EF_onlyMixtures_onlyGMA", {
##   inputs <- simulDat()
##   truth <- inputs$truth
##   dat <- inputs$dat
##   X <- inputs$X
##   Z_GMA <- inputs$Z_GMA

##   expVarComps <- c(
##     "GMA" = truth$var_GMA,
##     "Residual" = truth$var_error
##   )

##   ## fit the model
##   ## obsFit <- lmerGMASMA_EF(Response=dat$pheno, X=X, ListZ=list("GMA"=Z_GMA),
##   ##                         REML=TRUE)
##   ## if(FALSE){
##   ##   plot(x=lme4::ranef(obsFit)$GMA[,1], y=truth$GMAs, las=1,
##   ##        xlab="BLUP(GMA)", ylab="GMA")
##   ##   abline(a=0, b=1)
##   ## }
##   ## obsVarComps <- as.data.frame(lme4::VarCorr(obsFit))
##   ## obsVarComps <- c("GMA"=obsVarComps$vcov[1],
##   ##                  "Residual"=obsVarComps$vcov[2])
##   ## obsGMAs <- ranefGMASMA_EF(obsFit) # returns an error

##   ## check the results -> errors
##   ## expect_equal(obsVarComps["GMA"], expVarComps["GMA"], tolerance=2)
##   ## expect_equal(obsVarComps["Residual"], expVarComps["Residual"], tolerance=0.2)
##   ## expect_equal(cor(obsGMAs, truth$GMAs), 1, tolerance=0.1)
## })

test_that("lmerZ_onlyMixtures_onlyGMA", {
  inputs <- simulDat()
  truth <- inputs$truth
  dat <- inputs$dat
  Z_GMA <- inputs$Z_GMA

  expVarComps <- c(
    "GMA" = truth$var_GMA,
    "Residual" = truth$var_error
  )

  ## fit the model
  dat2 <- cbind(dat, "GMA" = dat$stands)
  myformula <- pheno ~ 1 + block + (1 | GMA)
  myblist <- list(list(
    "ff" = factor(colnames(Z_GMA)),
    "sm" = Matrix::Matrix(t(Z_GMA), sparse = TRUE),
    "nl" = as.integer(ncol(Z_GMA)),
    "cnms" = "(Intercept)"
  ))
  names(myblist) <- lme4:::barnames(lme4::findbars(lme4:::RHSForm(myformula)))
  obsFit <- lmerZ(
    formula = myformula, data = dat2,
    REML = TRUE, myblist = myblist
  )
  if (FALSE) {
    plot(
      x = lme4::ranef(obsFit)$GMA[, 1], y = truth$GMAs, las = 1,
      xlab = "BLUP(GMA)", ylab = "GMA"
    )
    abline(a = 0, b = 1)
  }
  obsVarComps <- as.data.frame(lme4::VarCorr(obsFit))
  obsVarComps <- c(
    "GMA" = obsVarComps$vcov[1],
    "Residual" = obsVarComps$vcov[2]
  )
  obsGMAs <- lme4::ranef(obsFit)$GMA[, 1]

  ## check the results
  expect_equal(obsVarComps["GMA"], expVarComps["GMA"], tolerance = 2)
  expect_equal(obsVarComps["Residual"], expVarComps["Residual"], tolerance = 0.2)
  expect_equal(cor(obsGMAs, truth$GMAs), 1, tolerance = 0.1)
})

test_that("lmerGMASMA_onlyMixtures_onlyGMA", {
  inputs <- simulDat()
  truth <- inputs$truth
  dat <- inputs$dat
  Z_GMA <- inputs$Z_GMA
  listContr <- inputs$listContr

  expVarComps <- c(
    "GMA" = truth$var_GMA,
    "Residual" = truth$var_error
  )

  ## fit the model
  myformFix <- pheno ~ 1 + block
  obsFit <- lmerGMASMA(myformFix, dat, list("GMA" = Z_GMA), listContr)
  if (FALSE) {
    plot(
      x = lme4::ranef(obsFit)$GMA[, 1], y = truth$GMAs, las = 1,
      xlab = "BLUP(GMA)", ylab = "GMA"
    )
    abline(a = 0, b = 1)
  }
  obsVarComps <- as.data.frame(lme4::VarCorr(obsFit))
  obsVarComps <- c(
    "GMA" = obsVarComps$vcov[1],
    "Residual" = obsVarComps$vcov[2]
  )
  obsGMAs <- lme4::ranef(obsFit)$GMA[, 1]

  ## check the results
  expect_equal(obsVarComps["GMA"], expVarComps["GMA"], tolerance = 2)
  expect_equal(obsVarComps["Residual"], expVarComps["Residual"], tolerance = 0.2)
  expect_equal(cor(obsGMAs, truth$GMAs), 1, tolerance = 0.1)
})
