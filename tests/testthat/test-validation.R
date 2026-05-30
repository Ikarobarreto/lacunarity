test_that("lac and genlac run on a valid binary series", {
  set.seed(1)
  x <- rbinom(500, 1, 0.7)
  expect_named(lac(x), c("y", "Ds", "s"))
  expect_named(genlac(x), c("s", "q", "yq", "Dqs"))
})

test_that("lac rejects invalid input with a clear message", {
  expect_error(lac(rnorm(100)), "0's and 1's")
  expect_error(lac(c(rbinom(99, 1, 0.5), NA)), "missing")
  expect_error(lac(matrix(rbinom(100, 1, 0.5), ncol = 2)), "vector")
  expect_error(lac(c(0, 1, 1)), "too short")
})
