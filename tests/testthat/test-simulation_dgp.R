context("Testing ATbounds: simulation_dgp")

test_that("ATE-oracle and the mean difference estimator should be similar when P(D=1|X) = 0.5", {
  set.seed(1)
  data <- simulation_dgp(10000, ps_spec = "overlap")
  y <- data$outcome
  d <- data$treat
  ate <- data$ate_oracle
  mde <- mean(d*y)/mean(d) - mean((1-d)*y)/mean(1-d)
  diff <- (abs(ate - mde) <= 1e-2)

  expect_true(diff)
  
})

test_that("ATE-oracle and the mean difference estimator should be different when P(D=1|X) is a function of X", {
  set.seed(1)
  data <- simulation_dgp(10000, ps_spec = "non-overlap")
  y <- data$outcome
  d <- data$treat
  ate <- data$ate_oracle
  mde <- mean(d*y)/mean(d) - mean((1-d)*y)/mean(1-d)
  diff <- (abs(ate - mde) <= 1e-2)
  
  expect_false(diff)
  
})

test_that("Discrete and continuous cases are different", {
  set.seed(1)
  data_d <- simulation_dgp(100, ps_spec = "non-overlap", x_discrete = TRUE)
  set.seed(1)
  data_c <- simulation_dgp(100, ps_spec = "non-overlap")
  expect_false(data_c$covariate[1] == data_d$covariate[1])
  
})