library(plantmix)
context("incomplete balanced design")

test_that("allocGenosMixsInterMixDesign_only-two-testers", {
  levGenos <- list(
    "geno_S1" = paste0("wheat", 1:4),
    "geno_S2" = c("pea1", "pea2")
  )
  dat_bl1 <- plantmix:::makeCompleteDesignOneBlockInterMixDesign(levGenos)
  levBlocks <- c("A", "B")
  datW <- plantmix:::makeCompleteDesignInterMixDesign(dat_bl1, levBlocks)
  datW <- datW[, c("geno_S1", "geno_S2", "type2", "block")]
  expected <- data.frame(
    geno_S1 = c(
      "wheat2", "wheat3", NA, NA, "wheat3", "wheat4",
      "wheat1", "wheat4", NA, NA, "wheat2", "wheat1"
    ),
    geno_S2 = c(
      NA, NA, "pea1", "pea2", "pea1", "pea2",
      NA, NA, "pea1", "pea2", "pea1", "pea2"
    ),
    type2 = rep(
      c("sole_S1", "sole_S1", "sole_S2", "sole_S2", "mix", "mix"),
      2
    ),
    block = rep(levBlocks, each = 6),
    stringsAsFactors = TRUE
  )
  expected$type2 <- factor(as.character(expected$type2),
    levels = c("sole_S1", "sole_S2", "mix")
  )
  observed <- plantmix:::allocGenosMixsInterMixDesign_OLD(
    datW, 2, levGenos,
    seed = 12345
  )
  idx <- order(observed$block, observed$type2, observed$geno_S2, observed$geno_S1)
  observed <- observed[idx, ]
  rownames(observed) <- NULL
  expect_equal(observed, expected)
  if (FALSE) {
    plantmix:::plotAllocSchemeInterMixDesign(observed)
  }

  observed2 <- plantmix:::allocGenosMixsInterMixDesign(
    datW, 2, levGenos,
    seed = 12345
  )
  idx <- order(observed2$block, observed2$type2, observed2$geno_S2, observed2$geno_S1)
  observed2 <- observed2[idx, ]
  rownames(observed2) <- NULL

  tmp <- table(observed2$block)
  expect_true(length(unique(tmp)) == 1)
  tmp <- table(observed2$block, observed2$type2)
  expect_true(all(apply(tmp, 2, function(x) {
    length(unique(x))
  }) == 1))
  for (levG2 in levels(datW$geno_S2)) {
    idx <- which(observed2$geno_S2 == levG2)
    tmp <- table(observed2$block[idx], observed2$type2[idx])
    expect_true(all(tmp[, "sole_S1"] == 0))
    expect_true(all(apply(tmp, 2, function(x) {
      length(unique(x))
    }) == 1))
  }
  for (levG1 in levels(datW$geno_S1)) {
    idx <- which(observed2$geno_S1 == levG1)
    tmp <- table(observed2$block[idx], observed2$type2[idx])
    expect_true(all(tmp[, "sole_S2"] == 0))
    tmp <- tmp[, -which(colnames(tmp) == "sole_S2")]
    expect_true(all(rowSums(tmp) == 1))
  }
})

test_that("allocGenosMixsInterMixDesign_more-than-two-testers", {
  levGenos <- list(
    "geno_S1" = paste0("wheat", 1:6),
    "geno_S2" = c("pea1", "pea2", "pea3")
  )
  dat_bl1 <- plantmix:::makeCompleteDesignOneBlockInterMixDesign(levGenos)
  levBlocks <- c("A", "B")
  datW <- plantmix:::makeCompleteDesignInterMixDesign(dat_bl1, levBlocks)
  datW <- datW[, c("geno_S1", "geno_S2", "type2", "block")]

  observed2 <- plantmix:::allocGenosMixsInterMixDesign(
    datW, 2, levGenos,
    seed = 12345
  )
  idx <- order(observed2$block, observed2$type2, observed2$geno_S2, observed2$geno_S1)
  observed2 <- observed2[idx, ]
  rownames(observed2) <- NULL

  tmp <- table(observed2$block)
  expect_true(length(unique(tmp)) == 1)
  tmp <- table(observed2$block, observed2$type2)
  expect_true(all(apply(tmp, 2, function(x) {
    length(unique(x))
  }) == 1))
  for (levG2 in levels(datW$geno_S2)) {
    idx <- which(observed2$geno_S2 == levG2)
    tmp <- table(observed2$block[idx], observed2$type2[idx])
    expect_true(all(tmp[, "sole_S1"] == 0))
    expect_true(all(apply(tmp, 2, function(x) {
      length(unique(x))
    }) == 1))
  }
  for (levG1 in levels(datW$geno_S1)) {
    idx <- which(observed2$geno_S1 == levG1)
    tmp <- table(observed2$block[idx], observed2$type2[idx])
    expect_true(all(tmp[, "sole_S2"] == 0))
    tmp <- tmp[, -which(colnames(tmp) == "sole_S2")]
    expect_true(all(rowSums(tmp) == 1))
  }
  ## expect_equal(observed2, expected)
  if (FALSE) {
    plantmix:::plotAllocSchemeInterMixDesign(observed2)
  }
})
