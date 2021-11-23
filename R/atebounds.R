#' @title Bounding the average treatment effect (ATE) 
#'
#' @description Bounds the average treatment effect (ATE) under the unconfoundedness assumption without the overlap condition.

#' @param Y n-dimensional vector of binary outcomes
#' @param D n-dimensional vector of binary treatments
#' @param X n by p matrix of covariates
#' @param rps n-dimensional vector of the reference propensity score
#' @param Q bandwidth parameter that determines the maximum number of observations for pooling information (default: Q = 3)
#' @param studentize TRUE if the columns of X are studentized and FALSE if not (default: TRUE)
#' @param alpha (1-alpha) nominal coverage probability for the confidence interval of ATE (default: 0.05)
#' @param x_discrete TRUE if the distribution of X is discrete and FALSE otherwise (default: FALSE)
#' @param n_hc number of hierarchical clusters to discretize non-discrete covariates; relevant only if x_discrete is FALSE.
#' The default choice is n_hc = ceiling(length(Y)/10), so that there are 10 observations in each cluster on average. 
#' 
#' @return An S3 object of type "ATbounds". The object has the following elements.
#' \item{call}{a call in which all of the specified arguments are specified by their full names}
#' \item{type}{ATE}
#' \item{cov_prob}{Confidence level: 1-alpha}
#' \item{y1_lb}{estimate of the lower bound on the average of Y(1), i.e. E[Y(1)]}
#' \item{y1_ub}{estimate of the upper bound on the average of Y(1), i.e. E[Y(1)]}
#' \item{y0_lb}{estimate of the lower bound on the average of Y(0), i.e. E[Y(0)]}
#' \item{y0_ub}{estimate of the upper bound on the average of Y(0), i.e. E[Y(0)]}
#' \item{est_lb}{estimate of the lower bound on ATE, i.e. E[Y(1) - Y(0)]}
#' \item{est_ub}{estimate of the upper bound on ATE, i.e. E[Y(1) - Y(0)]}
#' \item{est_rps}{the point estimate of ATE using the reference propensity score}
#' \item{se_lb}{standard error for the estimate of the lower bound on ATE}
#' \item{se_ub}{standard error for the estimate of the upper bound on ATE}
#' \item{ci_lb}{the lower end point of the confidence interval for ATE}
#' \item{ci_ub}{the upper end point of the confidence interval for ATE}
#' 
#' @examples
#'   Y <- RHC[,"survival"]
#'   D <- RHC[,"RHC"]
#'   X <- RHC[,c("age","edu")]
#'   rps <- rep(mean(D),length(D))
#'   results_ate <- atebounds(Y, D, X, rps, Q = 3)
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
atebounds <- function(Y, D, X, rps, Q = 3L, studentize = TRUE, alpha = 0.05, x_discrete = FALSE, n_hc = NULL){

  call <- match.call()
  
  X <- as.matrix(X)
  if (studentize == TRUE){
    X <- scale(X)  # centers and scales the columns of X
  }  

  n <- nrow(X)
  ymin <- min(Y)
  ymax <- max(Y)
  
  if (is.null(n_hc) == TRUE){ 
    n_hc = ceiling(n/10)  # number of clusters  
  }

  ### ATE estimation using reference propensity scores ###  
      
  y1_rps <- mean(D*Y/rps) 
  y0_rps <- mean((1-D)*Y/(1-rps))
  ate_rps <- y1_rps - y0_rps

  if (x_discrete == FALSE){  # Computing weights with non-discrete covariates

    hc <- stats::hclust(stats::dist(X), method = "complete")  # hierarchical cluster
    ind_Xunique <- stats::cutree(hc, k = n_hc)                # An index vector that has the same dimension as that of X
    mx <- n_hc
    
  } else if (x_discrete == TRUE){  # Computing weights with discrete covariates

    Xunique <- mgcv::uniquecombs(X)       # A matrix of unique rows from X
    ind_Xunique <- attr(Xunique,"index")  # An index vector that has the same dimension as that of X
    mx <- nrow(Xunique)                   # number of unique rows 

  } else {
    stop("x_discrete must be either TRUE or FALSE.") 
  }
  
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
      
      # Computing weights

      if (Q >= 1){
      
          qq <- min(Q,nx)
          k_upper <- 2*floor(qq/2)
          v_x1 <- 1
          v_x0 <- 1
          
          rps_x <- rps[disc_ind]
          
          if (x_discrete == FALSE){ 
            rps_x <- mean(rps_x)  # if X is non-discrete, take the average
          } else if (x_discrete == TRUE){
              rps_x <- unique(rps_x)
              if (length(rps_x) > 1){
                stop("The reference propensity score should be unique for the same value of X if X is discrete.")   
              }
          }    
          
          for (k in 0:k_upper){
              
               px1k <- ((rps_x-1)/rps_x)^k
               px0k <- (rps_x/(rps_x-1))^k
              
               if ((qq %% 2) == 1){ # if min(q,nx) is odd
                
                 term_x1 <- ((nx - nx1)/nx)*(1/choose(nx-1,qq-1))*choose(nx1,k)*choose(nx-1-nx1,qq-1-k)
                 term_x0 <- ((nx - nx0)/nx)*(1/choose(nx-1,qq-1))*choose(nx0,k)*choose(nx-1-nx0,qq-1-k)
                
               } else if ((qq %% 2) == 0){ # if min(q,nx) is even
                
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
      
      } else{
        stop("'Q' must be a positive integer.")      
      }
      
      res[i,1] <- ((mx*nx)/n)*(ymin + v_x1*(y1bar - ymin))
      res[i,2] <- ((mx*nx)/n)*(ymin + v_x0*(y0bar - ymin))
      res[i,3] <- ((mx*nx)/n)*(ymax + v_x1*(y1bar - ymax))
      res[i,4] <- ((mx*nx)/n)*(ymax + v_x0*(y0bar - ymax))
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
  
  # Stoye (2020) construction
  two_sided <- 1-alpha/2
  cv_norm <- stats::qnorm(two_sided)
  ci1_lb <- ate_lb - cv_norm*se_lb
  ci1_ub <- ate_ub + cv_norm*se_ub
  
  ate_star <- (se_ub*ate_lb + se_lb*ate_ub)/(se_lb + se_ub)
  se_star <- 2*(se_lb*se_ub)/(se_lb + se_ub) # This corresponds to rho=1 in Stoye (2020)
  ci2_lb <- ate_star - cv_norm*se_star
  ci2_ub <- ate_star + cv_norm*se_star
  
  if (ci1_lb <= ci1_ub){ # if the first confidence interval is non-empty
    ci_lb <- min(ci1_lb,ci2_lb)
    ci_ub <- max(ci1_ub,ci2_ub)
  } else {
    ci_lb <- ci2_lb
    ci_ub <- ci2_ub  
  }
  
  outputs = list(call = call, type = "ATE", cov_prob = (1-alpha),
                 "y1_lb"=y1_lb,"y1_ub"=y1_ub,"y0_lb"=y0_lb,"y0_ub"=y0_ub,
                 "est_lb"=ate_lb,"est_ub"=ate_ub,"est_rps"=ate_rps,
                 "se_lb"=se_lb,"se_ub"=se_ub,"ci_lb"=ci_lb,"ci_ub"=ci_ub)  

  class(outputs) = 'ATbounds'

  outputs
}
