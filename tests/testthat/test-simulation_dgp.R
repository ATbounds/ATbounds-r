context("Testing ATbounds")

test_that("ATE-oracle and the mean difference estimator should be similar when P(D=1|X) = 0.5", {
  
  data <- simulation_dgp(10000, ps_spec = "overlap")
  y <- data$outcome
  d <- data$treat
  ate <- data$ate_oracle
  mde <- mean(d*y)/mean(d) - mean((1-d)*y)/mean(1-d)
  diff <- (abs(ate - mde) <= 1e-2)

  expect_true(diff)
  
})

test_that("ATE-oracle and the mean difference estimator should be different when P(D=1|X) is a function of X", {
  
  data <- simulation_dgp(10000, ps_spec = "non-overlap")
  y <- data$outcome
  d <- data$treat
  ate <- data$ate_oracle
  mde <- mean(d*y)/mean(d) - mean((1-d)*y)/mean(1-d)
  diff <- (abs(ate - mde) <= 1e-2)
  
  expect_false(diff)
  
})