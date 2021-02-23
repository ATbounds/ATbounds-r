#' @title Bounding the average treatment effect on the treated (ATT) with discrete covariates 
#'
#' @description Bounds the average treatment effect on the treated (ATT) under the unconfoundedness assumption without the overlap condition.
#' This command is for the case when all the covariates are discrete.
#'
#' @param Y n-dimensional vector of binary outcomes
#' @param D n-dimensional vector of binary treatments
#' @param X n by p matrix of covariates
#' @param rps n-dimensional vector of the reference propensity score
#' @param Q bandwidth parameter that determines the maximum number of observations for pooling information (default: Q = 3)
#' @param studentize TRUE if X is studentized elementwise and FALSE if not (default: TRUE)
#' @param alpha (1-alpha) nominal coverage probability for the confidence interval of ATE (default: 0.05)
#' 
#' @return An S3 object of type "ATbounds". The object has the following elements.
#' \item{att_lb}{estimate of the lower bound on ATT, i.e. E[Y(1) - Y(0) | D = 1]}
#' \item{att_ub}{estimate of the upper bound on ATT, i.e. E[Y(1) - Y(0) | D = 1]}
#' \item{att_rps}{the point estimate of ATT using the reference propensity score}
#' \item{se_lb}{standard error for the estimate of the lower bound on ATT}
#' \item{se_ub}{standard error for the estimate of the upper bound on ATT}
#' \item{ci_lb}{the lower end point of the confidence interval for ATT}
#' \item{ci_ub}{the upper end point of the confidence interval for ATT}
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
attbounds_discrete_x <- function(Y, D, X, rps, Q = 3L, studentize = TRUE, alpha = 0.05){

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
  

  Xunique <- mgcv::uniquecombs(X)       # A matrix of unique rows from X
  ind_Xunique <- attr(Xunique,"index")  # An index vector that the same dimension as that of X
      
  mx <- nrow(Xunique) # number of unique rows 
      
  res <- matrix(NA,nrow=mx,ncol=2)
      
  for (i in 1:mx){ # this loop may not be fast if mx is very large
    
      disc_ind <- (ind_Xunique == i) 
      nx <- sum(disc_ind)    # number of obs. such that X_i = x for each row of x
      nx1 <- sum(D*disc_ind) # number of obs. such that X_i = x and D_i = 1 for each row of x
      nx0 <- nx - nx1        # number of obs. such that X_i = x and D_i = 0 for each row of x
      
      nx1 <- nx1 + (nx1 == 0) # replace nx1 with 1 when it is zero to avoid NaN 
      nx0 <- nx0 + (nx0 == 0) # replace nx0 with 1 when it is zero to avoid NaN 
      
      y1bar <- (sum(D*Y*disc_ind)/nx1)     # Dividing by zero never occurs because the numerator is zero whenever nx1 is zero  
      y0bar <- (sum((1-D)*Y*disc_ind)/nx0) # Dividing by zero never occurs because the numerator is zero whenever nx0 is zero 
      
      
      # Computing weights

      if (Q >= 1){

          qq <- min(Q,nx)
          k_upper <- 2*floor(qq/2)
          v_x0 <- nx1/nx
          
          for (k in 0:k_upper){

              rps_x <- rps[disc_ind]
              rps_x <- unique(rps_x)
              
              if (length(rps_x) > 1){
                stop("The reference propensity score should be unique for the same value of x.")   
              }          
              
              px0k <- (rps_x/(rps_x-1))^k
              #omega0 <- px0k
              
              if ((qq %% 2) == 1){ # if min(q,nx) is odd
                
                term_x0 <- ((nx - nx0)/nx)*(1/choose(nx-1,qq-1))*choose(nx0,k)*choose(nx-1-nx0,qq-1-k)
                
              } else if ((qq %% 2) == 0){ # min(q,nx) is even
                
                term_x0 <- (1/choose(nx,qq))*choose(nx0,k)*choose(nx-nx0,qq-k)
                
              } else{
                stop("'min(q,nx)' must be a positive integer.")      
              }  
              
              omega0 <- px0k*term_x0

              v_x0 <- v_x0 - omega0
            
          }

       } else{
         stop("'Q' must be a positive integer.")      
       }
        
       res[i,1] <- ((mx*nx)/n)*((nx1/nx)*(y1bar - ymax) - v_x0*(y0bar - ymax))
       res[i,2] <- ((mx*nx)/n)*((nx1/nx)*(y1bar - ymin) - v_x0*(y0bar - ymin))
        
   } 
    
  
  ### Obtain bound estimates ###
  
  est <- apply(res,2,mean)
  att_lb <- est[1]/mean(D)
  att_ub <- est[2]/mean(D)

  # Standard errors, while treating mean(D) fixed
  se <- apply(res,2,stats::sd)
  se_lb <- se[1]/(sqrt(mx)*mean(D))
  se_ub <- se[2]/(sqrt(mx)*mean(D))

  # Stoye (2020) construction
  two_sided <- 1-alpha/2
  cv_norm <- stats::qnorm(two_sided)
  ci1_lb <- att_lb - cv_norm*se_lb
  ci1_ub <- att_ub + cv_norm*se_ub

  att_star <- (se_ub*att_lb + se_lb*att_ub)/(se_lb + se_ub)
  se_star <- 2*(se_lb*se_ub)/(se_lb + se_ub) # This corresponds to rho=1 in Stoye (2020)
  ci2_lb <- att_star - cv_norm*se_star
  ci2_ub <- att_star + cv_norm*se_star

  if (ci1_lb <= ci1_ub){ # if the first confidence interval is non-empty
    ci_lb <- min(ci1_lb,ci2_lb)
    ci_ub <- max(ci1_ub,ci2_ub)
  } else {
    ci_lb <- ci2_lb
    ci_ub <- ci2_ub  
  }
  
  outputs = list("att_lb"=att_lb,"att_ub"=att_ub,"att_rps"=att_rps,
                 "se_lb"=se_lb,"se_ub"=se_ub,"ci_lb"=ci_lb,"ci_ub"=ci_ub)  


  class(outputs) = 'ATbounds'

  outputs
}

