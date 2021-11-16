context("Testing ATbounds::atebounds")

Y <- RHC[,"survival"]
D <- RHC[,"RHC"]
age <- round(RHC[,"age"])
female <- RHC[,"sex_Female"]
X <- cbind(age,female)
rps <- rep(mean(D),length(D))

test_that("Q must be strictly positive", {
  
  expect_error(atebounds(Y, D, X, rps, Q = 0))
  
})

test_that("Checking the summary results", {
  
  sum = summary(atebounds(Y, D, X, rps))
  expect_named(sum, c("Lower_Bound","Upper_Bound"))
  
})

test_that("'x_discrete' must be either TRUE or FALSE", {
  
  expect_error(atebounds(Y, D, X, rps, x_discrete=Maybe))
  
})

Y_EFM <- EFM[,"cesarean"]
D_EFM <- EFM[,"monitor"]
X_EFM <- as.matrix(EFM[,c("arrest", "breech", "nullipar", "year")])
rps_EFM <- rep(mean(D_EFM),length(D_EFM))  

test_that("Different Q generates different results", {
  
  res1 <- atebounds(Y_EFM, D_EFM, X_EFM, rps_EFM, Q=3, x_discrete=TRUE)
  res2 <- atebounds(Y_EFM, D_EFM, X_EFM, rps_EFM, Q=4, x_discrete=TRUE)
  
  expect_false(res1$ci_lb == res2$ci_lb)
  
})