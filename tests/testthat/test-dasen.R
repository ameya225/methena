test_that("memory-efficient dasen reproduces wateRmelon::dasen", {
  skip_if_not_installed("wateRmelon"); skip_if_not_installed("preprocessCore")
  set.seed(3); p <- 3000; n <- 12
  M <- matrix(exp(rnorm(p * n, 8, 1)), p, n); U <- matrix(exp(rnorm(p * n, 8, 1)), p, n)
  ids <- paste0("200000000", sprintf("%02d", 1:n), "_R0", rep(1:6, 2), "C0", rep(1:2, each = 6))
  colnames(M) <- colnames(U) <- ids; rownames(M) <- rownames(U) <- paste0("cg", 1:p)
  dv <- sample(1:2, p, TRUE, prob = c(0.3, 0.7))
  ref <- suppressWarnings(wateRmelon::dasen(M, U, onetwo = ifelse(dv == 1, "I", "II")))
  M2 <- methena:::dasen_channel(M, dv, methena:::sentrix_roco(ids))
  U2 <- methena:::dasen_channel(U, dv, NULL)
  mine <- M2 / (M2 + U2 + 100)
  expect_lt(max(abs(mine - ref)), 1e-8)
})
