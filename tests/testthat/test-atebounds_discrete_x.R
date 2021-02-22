context("Testing ATbounds")

test_that("Q must be strictly positive", {
  
  Y <- RHC[,"survival"]
  D <- RHC[,"RHC"]
  age <- round(RHC[,"age"])
  female <- RHC[,"sex_Female"]
  X <- cbind(age,female)
  rps <- rep(mean(D),length(D))

  expect_error(atebounds_discrete_x(Y, D, X, rps, Q = 0))
  
})