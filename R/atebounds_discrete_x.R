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
#' 
#' @return An S3 object of type "ATbounds". The object has the following elements.
#' \item{y1_lb}{the lower bound on the average of Y(1), i.e. E[Y(1)]}
#' \item{y1_ub}{the upper bound on the average of Y(1), i.e. E[Y(1)]}
#' \item{y0_lb}{the lower bound on the average of Y(0), i.e. E[Y(0)]}
#' \item{y0_ub}{the upper bound on the average of Y(0), i.e. E[Y(0)]}
#' \item{ate_lb}{estimate of the lower bound on ATE, i.e. E[Y(1) - Y(0)]}
#' \item{ate_ub}{estimate of the upper bound on ATE, i.e. E[Y(1) - Y(0)]}
#' \item{ate_rps}{the point estimate of ATE using the reference propensity scores}
#' \item{se_lb}{standard error for the estimate of the lower bound on ATE, i.e. E[Y(1) - Y(0)]}
#' \item{se_ub}{standard error for the estimate of the upper bound on ATE, i.e. E[Y(1) - Y(0)]}
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
atebounds_discrete_x <- function(Y, D, X, rps, Q = 2L, studentize = TRUE){

  X <- as.matrix(X)
  n <- nrow(X)
  
  ymin <- min(Y)
  ymax <- max(Y)

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
      
      Xunique <- mgcv::uniquecombs(X)      # A matrix of unique rows from X
      ind_Xunique <- attr(Xunique,"index")  # An index vector that the same dimension as that of X
      
      mx <- nrow(Xunique) # number of unique rows 
      
      res <- matrix(NA,nrow=mx,ncol=4)
      
      for (i in 1:mx){ # this loop may not be fast if mx is very large
        
        disc_ind <- (ind_Xunique == i) 
        nx <- sum(disc_ind)    # number of obs. such that X_i = x for each row of x
        nx1 <- sum(D*disc_ind) # number of obs. such that X_i = x and D_i = 1 for each row of x
        nx0 <- nx - nx1        # number of obs. such that X_i = x and D_i = 0 for each row of x
        
        nx1 <- nx1 + (nx1 == 0) # replace nx1 with 1 when it is zero to avoid NaN 
        nx0 <- nx0 + (nx0 == 0) # replace nx0 with 1 when it is zero to avoid NaN 
        
        y1bar <- (sum(D*Y*disc_ind)/nx1)     # Dividing by zero never occurs because the numerator is zero whenever nx1 is zero  
        y0bar <- (sum((1-D)*Y*disc_ind)/nx0) # Dividing by zero never occurs because the numerator is zero whenever nx0 is zero 
        
        qq <- min(Q,nx)
        k_upper <- 2*floor(qq/2)
        v_x1 <- 1
        v_x0 <- 1
        
        for (k in 0:k_upper){
          
           rps_x <- rps[disc_ind]
           rps_x <- unique(rps_x)
          
           if (length(rps_x) > 1){
             stop("The reference propensity score should be unique for the same value of x.")   
           }
        
          px1k <- ((rps_x-1)/rps_x)^k
          px0k <- (rps_x/(rps_x-1))^k
          
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
        
        res[i,1] <- ((mx*nx)/n)*(ymin + v_x1*(y1bar - ymin))
        res[i,2] <- ((mx*nx)/n)*(ymin + v_x0*(y0bar - ymin))
        res[i,3] <- ((mx*nx)/n)*(ymax + v_x1*(y1bar - ymax))
        res[i,4] <- ((mx*nx)/n)*(ymax + v_x0*(y0bar - ymax))
      } 
      
    } else{
      stop("'Q' must be a positive integer.")      
    }

  ### Obtain bound estimates ###
  
  est <- apply(res,2,mean)
  y1_lb <- est[1]
  y0_lb <- est[2]
  y1_ub <- est[3]
  y0_ub <- est[4]
  
  Lx <- res[,1]-res[,4]
  Ux <- res[,3]-res[,2]
  ate_lb <-  mean(Lx)
  ate_ub <-  mean(Ux)
    
  se_lb <- stats::sd(Lx)/sqrt(mx)
  se_ub <- stats::sd(Ux)/sqrt(mx)
  
  outputs = list("y1_lb"=y1_lb,"y1_ub"=y1_ub,"y0_lb"=y0_lb,"y0_ub"=y0_ub,
                 "ate_lb"=ate_lb,"ate_ub"=ate_ub,"ate_rps"=ate_rps,"se_lb"=se_lb,"se_ub"=se_ub)  

  class(outputs) = 'ATbounds'

  outputs
}

