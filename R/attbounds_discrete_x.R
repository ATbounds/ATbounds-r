#' @title Bounding the average treatment effect on the treated (ATT) with discrete covariates 
#'
#' @description Bounds the average treatment effect on the treated (ATT) under the unconfoundedness assumption without the overlap condition.
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
#' \item{att_lb}{the lower bound on ATT, i.e. E[Y(1) - Y(0) | D = 1]}
#' \item{att_ub}{the upper bound on ATT, i.e. E[Y(1) - Y(0) | D = 1]}
#' \item{att_rps}{the point estimate of ATT using the reference propensity scores}
#' 
#' @examples
#'   Y <- RHC[,"survival"]
#'   D <- RHC[,"RHC"]
#'   age <- round(RHC[,"age"])
#'   female <- RHC[,"sex_Female"]
#'   X <- cbind(age,female)
#'   rps <- rep(mean(D),length(D))
#'   results_att <- attbounds_discrete_x(Y, D, X, rps, Q = 3)
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
attbounds_discrete_x <- function(Y, D, X, rps, Q = 2L, studentize = TRUE){

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
  
  ### ATT estimation using reference propensity scores  ###    
  
  rps_wt <- rps/(1-rps)      
  att_rps <- sum(D*Y-rps_wt*(1-D)*Y)/sum(D)      

  ### Computing weights with discrete covariates ###
  
    if (Q == 1){

      att_wt <- 0

    } else if (Q > 1){
      
      Xunique <- mgcv::uniquecombs(X)      # A matrix of unique rows from X
      ind_Xunique <- attr(Xunique,"index")  # An index vector that the same dimension as that of X
      
      mx <- nrow(Xunique) # number of unique rows 
      
      res <- matrix(NA,nrow=mx,ncol=2)
      
      for (i in 1:mx){ # this loop may not be fast if mx is very large
        
        disc_ind <- (ind_Xunique == i) 
        nx <- sum(disc_ind)    # number of obs. such that X_i = x for each row of x
        nx1 <- sum(D*disc_ind) # number of obs. such that X_i = x and D_i = 1 for each row of x
        nx0 <- nx - nx1        # number of obs. such that X_i = x and D_i = 0 for each row of x
        
        nx0 <- nx0 + (nx0 == 0) # replace nx0 with 1 when it is zero to avoid NaN 
        
        y0bar <- (sum((1-D)*Y*disc_ind)/nx0) # Dividing by zero never occurs because the numerator is zero whenever nx0 is zero 
        
        qq <- min(Q,nx)
        k_upper <- 2*floor(qq/2)
        v_x0 <- nx1
        
        for (k in 0:k_upper){

          rps_x <- rps[disc_ind]
          rps_x <- unique(rps_x)
          
          if (length(rps_x) > 1){
            stop("The reference propensity score should be unique for the same value of x.")   
          }          
          
          px0k <- (rps_x/(rps_x-1))^k
          omega0 <- px0k
          
          if ((qq %% 2) == 1){ # if min(q,nx) is odd
            
            term_x0 <- ((nx - nx0)/nx)*(1/choose(nx-1,qq-1))*choose(nx0,k)*choose(nx-1-nx0,qq-1-k)
            
          } else if ((qq %% 2) == 0){ # min(q,nx) is even
            
            term_x0 <- (1/choose(nx,qq))*choose(nx0,k)*choose(nx-nx0,qq-k)
            
          } else{
            stop("'min(q,nx)' must be a positive integer.")      
          }  
          
          omega0 <- px0k*term_x0

          v_x0 <- v_x0 - nx*omega0
          
        }
      
        res[i,1] <- (mx/n)*(v_x0*(y0bar - ymax))
        res[i,2] <- (mx/n)*(v_x0*(y0bar - ymin))
        
      } 
      
    } else{
      stop("'Q' must be a positive integer.")      
    }

  ### Obtain bound estimates ###
  
  est <- apply(res,2,mean)
  att_lb <- (sum(D*Y)/sum(D) - ymax) - est[1]/mean(D)
  att_ub <- (sum(D*Y)/sum(D) - ymin) - est[2]/mean(D)
  
  outputs = list("att_lb"=att_lb,"att_ub"=att_ub,"att_rps"=att_rps)

  class(outputs) = 'ATbounds'

  outputs
}

