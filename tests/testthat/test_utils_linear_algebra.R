library(plantmix)

test_that("vec_without-names", {
  mat <- matrix(
    c(
      1, 2, 3,
      4, 5, 6
    ),
    nrow = 2, ncol = 3, byrow = TRUE
  )
  expected <- c(1, 4, 2, 5, 3, 6)
  observed <- vec(mat)
  expect_equal(observed, expected)
})

test_that("vec_with-names", {
  mat <- matrix(
    c(
      1, 2, 3,
      4, 5, 6
    ),
    nrow = 2, ncol = 3, byrow = TRUE,
    dimnames = list(
      as.character(1:2),
      letters[1:3]
    )
  )
  expected <- setNames(
    c(1, 4, 2, 5, 3, 6),
    c(
      "1-a", "2-a",
      "1-b", "2-b",
      "1-c", "2-c"
    )
  )
  observed <- vec(mat)
  expect_equal(observed, expected)
})

test_that("invvec_without-names", {
  x <- c(1, 4, 2, 5, 3, 6)
  expected <- matrix(
    c(
      1, 2, 3,
      4, 5, 6
    ),
    nrow = 2, ncol = 3, byrow = TRUE
  )
  observed <- invvec(x, 3)
  expect_equal(observed, expected)
})

test_that("invvec_without-names", {
  x <- setNames(
    c(1, 4, 2, 5, 3, 6),
    c("1-a", "2-a", "1-b", "2-b", "1-c", "2-c")
  )
  expected <- matrix(
    c(
      1, 2, 3,
      4, 5, 6
    ),
    nrow = 2, ncol = 3, byrow = TRUE,
    dimnames = list(
      as.character(1:2),
      letters[1:3]
    )
  )
  observed <- invvec(x, 3, sep = "-")
  expect_equal(observed, expected)
})

test_that("invvecMixes", {
  datL <- data.frame(
    ID = c("g1+g2", "g1+g2", "g1+g2", "g1+g2"),
    focal = c("g1", "g2", "g1", "g2"),
    neighbor = c("g2", "g1", "g2", "g1"),
    block = c("A", "A", "B", "B"),
    yield = c(1, 2, 3, 4)
  )
  expected <- matrix(
    c(
      1, 2,
      3, 4
    ),
    nrow = 2, ncol = 2, byrow = TRUE
  )
  observed <- invvecMixes(datL$yield, 2)
  expect_equal(observed, expected)
})
