#' @title Bounding the average treatment effect on the treated (ATT)
#'
#' @description Bounds the average treatment effect on the treated (ATT) under the unconfoundedness assumption without the overlap condition
#'
#' @param Y n-dimensional vector of binary outcomes
#' @param D n-dimensional vector of binary treatments
#' @param X n by p matrix of covariates
#' @param rps n-dimensional vector of reference propensity scores
#' @param Q polynomial order that uses the (Q-1) nearest neighbors excluding own observations (default: Q = 2)
#' @param n_permute number of permutations to shuffle the data (default: 0)
#' @param studentize TRUE if x is studentized elementwise and FALSE if not (default: TRUE)
#' 
#' @return An S3 object of type "ATbounds". The object has the following elements.
#' \item{lb}{the lower bound on ATT, i.e. E[Y(1) - Y(0) | D = 1]}
#' \item{ub}{the upper bound on ATT, i.e. E[Y(1) - Y(0) | D = 1]}
#' \item{att_rps}{the point estimate of ATT using the reference propensity scores}
#' 
#' @examples
#'   Y <- RHC[,"survival"]
#'   D <- RHC[,"RHC"]
#'   X <- RHC[,-c(1,2)]
#'   rps <- rep(mean(D),length(D))  
#'   results_att <- attbounds(Y, D, X, rps, Q = 3)
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
simulation_dgp <- function(n, ps_spec = "overlap", X, rps, Q = 2L,x_discrete = TRUE, studentize = TRUE){
    
  x <- runif(n, min = -3, max = 3)
  
  if (x_discrete == TRUE){
      x <- round(x*10)/10 # discrete X in {-3.0, -2.9, ..., 3.0}  
  }
  
  if (ps_spec == "overlap"){  
    px <- 0.5
  } else if (ps_spec == "non-overlap"){
    px <- 0.25*(x >= 2) + 0.5*(abs(x) < 2) 
  }
  
  treat <- as.integer(px < runif(n)) # treat = 1 always if x <= -2 for the non-overlap case
  
  ps <- 1-px
  
  y_1 <- 1 + px + rnorm(n)
  y_0 <- px + rnorm(n)
  y_1  <- as.integer(y_1 > 0)
  y_0  <- as.integer(y_0 > 0)
  
  y <- treat*y_1 + (1-treat)*y_0
  
  y1 <- treat*y
  y0 <- (1-treat)*y
  
  ate_oracle <- mean(y1 - y0)
  att_oracle <- mean(treat*(y1-y0))/mean(ps)
  
  outputs <- list("outcome"=y,"treat"=treat,"covariate"=x,
                  "ate_oracle"=ate_oracle,"att_oracle"=att_oracle)
  
  class(outputs) = 'ATbounds'

outputs    
}
