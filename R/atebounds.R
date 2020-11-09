#' @title Bounding the average treatment effect
#'
#' @description Bounds the average treatment effect under the unconfoundedness assumption without the overlap condition
#'
#' @param y n-dimensional vector of binary outcomes
#' @param t n-dimensional vector of binary treatments
#' @param x n by p matrix of covariates
#' @param rps n-dimensional vector of reference propensity scores
#' @param q polynomial order (default: q = 2, which uses the nearest neighbor excluding own observations)
#' @param permute_max maximum number of permutations to shuffle the data (default: 0)
#' @param discrete TRUE if x includes only discrete covariates and FALSE if not (default: FALSE)
#' @param studentize TRUE if x is studentized elementwise and FALSE if not (default: TRUE)
#' @param small_c a small positive constant to determine the two covariate vectors are identical (default: 1e-8). 
#' This constatn is only used when the option 'discrete' is TRUE. 
#' 
#' @return An S3 object of type "ATbounds". The object has the following elements.
#' \item{y1_lb}{the lower bound of the average of Y(1), i.e. E[Y(1)]}
#' \item{y1_ub}{the upper bound of the average of Y(1), i.e. E[Y(1)]}
#' \item{y0_lb}{the lower bound of the average of Y(0), i.e. E[Y(0)]}
#' \item{y0_ub}{the upper bound of the average of Y(0), i.e. E[Y(0)]}
#' \item{ate_lb}{the lower bound of ATE, i.e. E[Y(1) - Y(0)]}
#' \item{ate_ub}{the upper bound of ATE, i.e. E[Y(1) - Y(0)]}
#' \item{ate_rps}{the point estimate of ATE using the reference propensity score}
#' 
#' @examples
#'   Y <- RHC[,"survival"]
#'   D <- RHC[,"RHC"]
#'   X <- RHC[,-c(1,2)]
#'   rps <- rep(mean(D),length(D))  
#'   results_ate <- atebounds(Y, D, X, rps, q = 3)
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
atebounds <- function(y, t, x, rps, q = 2L, permute_max = 0L, discrete = FALSE, studentize = TRUE, small_c = 1e-8){

  # Check whether "discrete" is TRUE or FALSE
  if ((discrete != TRUE) & (discrete != FALSE)){
    stop("'discrete' must be either TRUE or FALSE.")
  }
  
  # Studentize covariates elementwise
  
  x <- as.matrix(x)
  n <- nrow(x)
  
  if (studentize == TRUE){
  sd_x <- apply(x,2,stats::sd)
  sd_x <- matrix(sd_x,nrow=ncol(x),ncol=n)
  m_x <- apply(x,2,mean)
  m_x <- matrix(m_x,nrow=ncol(x),ncol=n)
  x <- (x-t(m_x))/t(sd_x) 
  }     
  
  results <- {}
  
  # Reorder the sample to break ties
  
  for (j in 0:permute_max){
  
    if ( permute_max > 0){    
    
      data <- cbind(y,t,x,rps)
      n <- nrow(data)  
      data <- data[sample.int(n,n),] # random permutation to shuffle the data  
      y <- data[,1]
      t <- data[,2]
      x <- data[,(3:(ncol(data)-1))]
      rps <- data[,ncol(data)]
    }
  
  ### ATE estimation using reference propensity scores ###  
      
  y1_rps <- mean(t*y/rps) 
  y0_rps <- mean((1-t)*y/(1-rps))
  ate_rps <- y1_rps - y0_rps

  ### Computing weights with discrete covariates ###
  
  if (discrete == TRUE){
    
    if (q == 1){
      y1_wt <- 1
      y0_wt <- 1
    } else if (q > 1){
      
      nn_data <- FNN::get.knnx(x, x, k=n) # this step may not be fast if n is large
      nn_i <- nn_data$nn.index
      nn_d <- nn_data$nn.dist # to get the Euclidiean distances among all covariates
      
      nn_t <- {} # n by q matrix of treatment values for counting nx1 and nx0
      for (k in 1:ncol(nn_i)){
        nn_t <- cbind(nn_t, t[nn_i[,k]])
      }
      
      nx <- rowSums(nn_d < small_c) # number of obs. such that X_i = x for each row of x
      nx1 <- rowSums(nn_t*(nn_d < small_c)) # number of obs. such that X_i = x and D_i = 1 for each row of x
      nx0 <- nx - nx1
      
      for (i in 1:n){
        
        qq <- min(q,nx[i])
        k_upper <- 2*floor(qq/2)
        v_x1 <- 0
        v_x0 <- 0
        
        for (k in 0:k_upper){
          
          px1k <- ((rps[i]-1)/rps[i])^k
          px0k <- (rps[i]/(rps[i]-1))^k
          
          omega1 <- px1k
          omega0 <- px0k
          
          if ((qq %% 2) == 1){ # if min(q,nx) is odd
            
            term_x1 <- ((nx[i] - nx1[i])/nx[i])*(1/choose(nx[i]-1,qq-1))*choose(nx1[i],k)*choose(nx[i]-1-nx1[i],qq-1-k)
            term_x0 <- ((nx[i] - nx0[i])/nx[i])*(1/choose(nx[i]-1,qq-1))*choose(nx0[i],k)*choose(nx[i]-1-nx0[i],qq-1-k)
            
          } else if ((q %% 2) == 0){ # min(q,nx) is even
            
            term_x1 <- (1/choose(nx[i],qq))*choose(nx1[i],k)*choose(nx[i]-nx1[i],qq-k)
            term_x0 <- (1/choose(nx[i],qq))*choose(nx0[i],k)*choose(nx[i]-nx0[i],qq-k)
            
          } else{
            stop("'q' must be a positive integer.")      
          }  
          
          omega1 <- px1k*term_x1
          omega0 <- px0k*term_x0

          v_x1 <- v_x1 + omega1
          v_x0 <- v_x0 + omega0
          
        }
        
      } 
      
      nx1c <- nx1 + (nx1 == 0L)
      nx0c <- nx0 + (nx0 == 0L)      
      y1_wt <- nx*v_x1/nx1c
      y0_wt <- nx*v_x0/nx0c
    } else{
      stop("'q' must be a positive integer.")      
    }
  }
  
  
  ### Computing weights with continuous covariates ###
  
  if (discrete == FALSE){
    
    if (q == 1){
      y1_wt <- 1
      y0_wt <- 1
    } else if (q > 1){
      
      nn_data <- FNN::get.knnx(x, x, k=q)
      nn_i <- nn_data$nn.index
      nn_d <- nn_data$nn.dist
      
      nn_t <- {} # n by q matrix of treatment values for NN estimation
      for (k in 1:ncol(nn_i)){
        nn_t <- cbind(nn_t, t[nn_i[,k]])
      }
      
      nx <- q
      nx1 <- rowSums(nn_t)
      nx0 <- nx - nx1
      
      px1 <- ((rps-1)/rps)^nx1
      px0 <- (rps/(rps-1))^nx0

      if ((q %% 2) == 1){ # if q is odd and q > 1
        v_x1 <- 1 - (nx0/nx)*px1
        v_x0 <- 1 - (nx1/nx)*px0
      } else if ((q %% 2) == 0){ # q is even
        v_x1 <- 1 - px1
        v_x0 <- 1 - px0
      } else{
        stop("'q' must be a positive integer.")      
      }  
      
      nx1c <- nx1 + (nx1 == 0L)
      nx0c <- nx0 + (nx0 == 0L)      
      y1_wt <- nx*v_x1/nx1c
      y0_wt <- nx*v_x0/nx0c
    } else{
      stop("'q' must be a positive integer.")      
    }
  }
  
  ### Obtain bound estimates ###
  
    y1_lb <- min(y) + y1_wt*(t == 1)*(y-min(y))
    y1_ub <- max(y) + y1_wt*(t == 1)*(y-max(y))
    y1_lb <- mean(y1_lb)
    y1_ub <- mean(y1_ub)
    
    y0_lb <- min(y) + y0_wt*(t == 0)*(y-min(y))
    y0_ub <- max(y) + y0_wt*(t == 0)*(y-max(y))
    y0_lb <- mean(y0_lb)
    y0_ub <- mean(y0_ub)
    
  result <- matrix(NA, nrow = 1, ncol = 4) 
  result <- c(y1_lb, y1_ub, y0_lb, y0_ub)
  results <- rbind(results,result)  
  } 

  
  est <- apply(results,2,mean)

  if (FALSE){
  y1_lb <- est[1]*(est[1] <= est[2]) + y1_rps*(est[1] > est[2])
  y1_ub <- est[2]*(est[1] <= est[2]) + y1_rps*(est[1] > est[2])
  y0_lb <- est[3]*(est[3] <= est[4]) + y0_rps*(est[3] > est[4])
  y0_ub <- est[4]*(est[3] <= est[4]) + y0_rps*(est[3] > est[4])
  }
  
  y1_lb <- est[1]
  y1_ub <- est[2]
  y0_lb <- est[3]
  y0_ub <- est[4]
  
  ate_lb <-  y1_lb - y0_ub
  ate_ub <-  y1_ub - y0_lb
  
  outputs = list("y1_lb"=y1_lb,"y1_ub"=y1_ub,
                 "y0_lb"=y0_lb,"y0_ub"=y0_ub,
                 "ate_lb"=ate_lb,"ate_ub"=ate_ub,
                 "ate_rps"=ate_rps)  

    class(outputs) = 'ATbounds'

    outputs
}

