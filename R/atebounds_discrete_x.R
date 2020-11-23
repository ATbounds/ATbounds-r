#' @title Bounding the average treatment effect (ATE) with discrete covariates 
#'
#' @description Bounds the average treatment effect (ATE) under the unconfoundedness assumption without the overlap condition.
#' This command is for the case when all the covariates are discrete.
#'
#' @param Y n-dimensional vector of binary outcomes
#' @param D n-dimensional vector of binary treatments
#' @param X n by p matrix of covariates
#' @param rps n-dimensional vector of reference propensity scores
#' @param Q polynomial order that uses the (Q-1) nearest neighbors excluding own observations (default: Q = 2)
#' @param studentize TRUE if X is studentized elementwise and FALSE if not (default: TRUE)
#' @param small_c a small positive constant to determine the two covariate vectors are identical (default: 1e-8)
#' 
#' @return An S3 object of type "ATbounds". The object has the following elements.
#' \item{y1_lb}{the lower bound on the average of Y(1), i.e. E[Y(1)]}
#' \item{y1_ub}{the upper bound on the average of Y(1), i.e. E[Y(1)]}
#' \item{y0_lb}{the lower bound on the average of Y(0), i.e. E[Y(0)]}
#' \item{y0_ub}{the upper bound on the average of Y(0), i.e. E[Y(0)]}
#' \item{ate_lb}{the lower bound on ATE, i.e. E[Y(1) - Y(0)]}
#' \item{ate_ub}{the upper bound on ATE, i.e. E[Y(1) - Y(0)]}
#' \item{ate_rps}{the point estimate of ATE using the reference propensity scores}
#' 
#' @examples
#'   Y <- RHC[,"survival"]
#'   D <- RHC[,"RHC"]
#'   age <- round(RHC[,"age"])
#'   female <- RHC[,"sex_Female"]
#'   X <- cbind(age,female)
#'   rps <- rep(mean(D),length(D))
#'   results_ate <- atebounds_discrete_x(Y, D, X, rps, Q = 3)
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
atebounds_discrete_x <- function(Y, D, X, rps, Q = 2L, studentize = TRUE, small_c = 1e-8){

  X <- as.matrix(X)
  n <- nrow(X)

  # Studentize covariates elementwise  
    
  if (studentize == TRUE){
  sd_x <- apply(X,2,stats::sd)
  sd_x <- matrix(sd_x,nrow=ncol(X),ncol=n)
  m_x <- apply(X,2,mean)
  m_x <- matrix(m_x,nrow=ncol(X),ncol=n)
  X <- (X-t(m_x))/t(sd_x) 
  }     
  
  ### ATE estimation using reference propensity scores ###  
      
  y1_rps <- mean(D*Y/rps) 
  y0_rps <- mean((1-D)*Y/(1-rps))
  ate_rps <- y1_rps - y0_rps

  ### Computing weights with discrete covariates ###
  
    if (Q == 1){
      y1_wt <- 1
      y0_wt <- 1
    } else if (Q > 1){
      
      y1_wt <- rep(NA,n)
      y0_wt <- rep(NA,n)
      
      for (i in 1:n){ # this loop may not be fast if n is very large
        
        xi <- t(matrix(X[i,],nrow=ncol(X),ncol=n))      
        dist <- sqrt(rowSums((X-xi)^2))
        dist_ind <- (dist < small_c)
        nx <- sum(dist_ind)    # number of obs. such that X_i = x for each row of x
        nx1 <- sum(D*dist_ind) # number of obs. such that X_i = x and D_i = 1 for each row of x
        nx0 <- nx - nx1        # number of obs. such that X_i = x and D_i = 0 for each row of x
        
        qq <- min(Q,nx)
        k_upper <- 2*floor(qq/2)
        v_x1 <- 1
        v_x0 <- 1
        
        for (k in 0:k_upper){
          
          px1k <- ((rps[i]-1)/rps[i])^k
          px0k <- (rps[i]/(rps[i]-1))^k
          
          omega1 <- px1k
          omega0 <- px0k
          
          if ((qq %% 2) == 1){ # if min(q,nx) is odd
            
            term_x1 <- ((nx - nx1)/nx)*(1/choose(nx-1,qq-1))*choose(nx1,k)*choose(nx-1-nx1,qq-1-k)
            term_x0 <- ((nx - nx0)/nx)*(1/choose(nx-1,qq-1))*choose(nx0,k)*choose(nx-1-nx0,qq-1-k)
            
          } else if ((qq %% 2) == 0){ # min(q,nx) is even
            
            term_x1 <- (1/choose(nx,qq))*choose(nx1,k)*choose(nx-nx1,qq-k)
            term_x0 <- (1/choose(nx,qq))*choose(nx0,k)*choose(nx-nx0,qq-k)
            
          } else{
            stop("'min(q,nx)' must be a positive integer.")      
          }  
          
          omega1 <- px1k*term_x1
          omega0 <- px0k*term_x0

          v_x1 <- v_x1 - omega1
          v_x0 <- v_x0 - omega0
          
        }
      
        # the second term (nx1 == 0L) is added to ensure that
        # y1_wt is well defined. However, nx1 is always 
        # strictly positive if D == 1. The same applies to 
        # (nx0 == 0L). 
        
        nx1c <- nx1 + (nx1 == 0L) 
        nx0c <- nx0 + (nx0 == 0L)      
        y1_wt[i] <- nx*v_x1/nx1c 
        y0_wt[i] <- nx*v_x0/nx0c  
      } 
      
    } else{
      stop("'Q' must be a positive integer.")      
    }

  ### Obtain bound estimates ###
  
    y1_lb <- min(Y) + y1_wt*(D == 1)*(Y-min(Y))
    y1_ub <- max(Y) + y1_wt*(D == 1)*(Y-max(Y))
    y1_lb <- mean(y1_lb)
    y1_ub <- mean(y1_ub)
    
    y0_lb <- min(Y) + y0_wt*(D == 0)*(Y-min(Y))
    y0_ub <- max(Y) + y0_wt*(D == 0)*(Y-max(Y))
    y0_lb <- mean(y0_lb)
    y0_ub <- mean(y0_ub)
    
    ate_lb <-  y1_lb - y0_ub
    ate_ub <-  y1_ub - y0_lb
  
    outputs = list("y1_lb"=y1_lb,"y1_ub"=y1_ub,"y0_lb"=y0_lb,"y0_ub"=y0_ub,
                   "ate_lb"=ate_lb,"ate_ub"=ate_ub,"ate_rps"=ate_rps)  

    class(outputs) = 'ATbounds'

    outputs
}

