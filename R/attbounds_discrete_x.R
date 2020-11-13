#' @title Bounding the average treatment effect on the treated with discrete covariates 
#'
#' @description Bounds the average treatment effect on the treated under the unconfoundedness assumption without the overlap condition.
#' This command is for the case when all the covariates are discrete.
#'
#' @param Y n-dimensional vector of binary outcomes
#' @param D n-dimensional vector of binary treatments
#' @param X n by p matrix of covariates
#' @param rps n-dimensional vector of reference propensity scores
#' @param Q maximal polynomial order that uses the (Q-1) nearest neighbors excluding own observations (default: Q = 2)
#' @param studentize TRUE if X is studentized elementwise and FALSE if not (default: TRUE)
#' @param small_c a small positive constant to determine the two covariate vectors are identical (default: 1e-8)
#' 
#' @return An S3 object of type "ATbounds". The object has the following elements.
#' \item{att_lb}{the lower bound of ATT, i.e. E[Y(1) - Y(0) | T = 1]}
#' \item{att_ub}{the upper bound of ATT, i.e. E[Y(1) - Y(0) | T = 1]}
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
attbounds_discrete_x <- function(Y, D, X, rps, Q = 2L, studentize = TRUE, small_c = 1e-8){

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
  
  ### ATT estimation using reference propensity scores  ###    
  
  rps_wt <- rps/(1-rps)      
  att_rps <- sum(D*Y-rps_wt*(1-D)*Y)/sum(D)      

  ### Computing weights with discrete covariates ###
  
    if (Q == 1){

      att_wt <- 0

    } else if (Q > 1){
      
      att_wt <- rep(NA,n)

      for (i in 1:n){ # this step may not be fast if n is very large
        
        xi <- t(matrix(X[i,],nrow=ncol(X),ncol=n))      
        dist <- sqrt(rowSums((X-xi)^2))
        dist_ind <- (dist < small_c)
        nx <- sum(dist_ind)    # number of obs. such that X_i = x for each row of x
        nx1 <- sum(D*dist_ind) # number of obs. such that X_i = x and D_i = 1 for each row of x
        nx0 <- nx - nx1        # number of obs. such that X_i = x and D_i = 0 for each row of x
        
        qq <- min(Q,nx)
        k_upper <- 2*floor(qq/2)
        v_x0 <- nx1
        
        for (k in 0:k_upper){
          
          px0k <- (rps[i]/(rps[i]-1))^k
          omega0 <- px0k
          
          if ((qq %% 2) == 1){ # if min(q,nx) is odd
            
            term_x0 <- ((nx - nx0)/nx)*(1/choose(nx-1,qq-1))*choose(nx0,k)*choose(nx-1-nx0,qq-1-k)
            
          } else if ((qq %% 2) == 0){ # min(q,nx) is even
            
            term_x0 <- (1/choose(nx,qq))*choose(nx0,k)*choose(nx-nx0,qq-k)
            
          } else{
            stop("'Q' must be a positive integer.")      
          }  
          
          omega0 <- px0k*term_x0

          v_x0 <- v_x0 - nx*omega0
          
        }
      
        nx0c <- nx0 + (nx0 == 0L)      

        att_wt[i] <- v_x0/nx0c  
        
      } 
      
    } else{
      stop("'Q' must be a positive integer.")      
    }

  ### Obtain bound estimates ###
  
  att_lb <- D*(Y-max(Y)) - att_wt*(1-D)*(Y-max(Y))
  att_ub <- D*(Y-min(Y)) - att_wt*(1-D)*(Y-min(Y))
  att_lb <- sum(att_lb)/sum(D)
  att_ub <- sum(att_ub)/sum(D)
  
  outputs = list("att_lb"=att_lb,"att_ub"=att_ub,"att_rps"=att_rps)

  class(outputs) = 'ATbounds'

  outputs
}

