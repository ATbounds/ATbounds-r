#' @title Bounding the average treatment effect (ATE)
#'
#' @description Bounds the average treatment effect (ATE) under the unconfoundedness assumption without the overlap condition.
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
#'   X <- RHC[,-c(1,2)]
#'   rps <- rep(mean(D),length(D))  
#'   results_ate <- atebounds(Y, D, X, rps, Q = 3)
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
atebounds <- function(Y, D, X, rps, Q = 2L, n_permute = 0L, studentize = TRUE){

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
  
  results <- {}
  
  # Reorder the sample to break ties
  
  for (j in 0:n_permute){
  
    if ( n_permute > 0){    
    
      data <- cbind(Y,D,X,rps)
      data <- data[sample.int(n,n),] # random permutation to shuffle the data  
      Y <- data[,1]
      D <- data[,2]
      X <- data[,(3:(ncol(data)-1))]
      rps <- data[,ncol(data)]
    }

  ### Computing weights with continuous and discrete covariates ###
  
    if (Q == 1){
      y1_wt <- 1
      y0_wt <- 1
    } else if (Q > 1){
      
      nn_data <- FNN::get.knnx(X, X, k=Q)
      nn_i <- nn_data$nn.index
      nn_d <- nn_data$nn.dist
      
      nn_t <- {} # n by Q matrix of treatment values for NN estimation
      for (k in 1:ncol(nn_i)){
        nn_t <- cbind(nn_t, D[nn_i[,k]])
      }
      
      nx <- Q
      nx1 <- rowSums(nn_t)
      nx0 <- nx - nx1
      
      px1 <- ((rps-1)/rps)^nx1
      px0 <- (rps/(rps-1))^nx0

      if ((Q %% 2) == 1){ # if Q is odd and Q > 1
        v_x1 <- 1 - (nx0/nx)*px1
        v_x0 <- 1 - (nx1/nx)*px0
      } else if ((Q %% 2) == 0){ # Q is even
        v_x1 <- 1 - px1
        v_x0 <- 1 - px0
      } else{
        stop("'Q' must be a positive integer.")      
      }  
      
      nx1c <- nx1 + (nx1 == 0L)
      nx0c <- nx0 + (nx0 == 0L)      
      y1_wt <- nx*v_x1/nx1c
      y0_wt <- nx*v_x0/nx0c
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
    
  result <- matrix(NA, nrow = 1, ncol = 4) 
  result <- c(y1_lb, y1_ub, y0_lb, y0_ub)
  results <- rbind(results,result)  
  } 

  est <- apply(results,2,mean)

  y1_lb <- est[1]
  y1_ub <- est[2]
  y0_lb <- est[3]
  y0_ub <- est[4]
  
  ate_lb <-  y1_lb - y0_ub
  ate_ub <-  y1_ub - y0_lb
  
  outputs = list("y1_lb"=y1_lb,"y1_ub"=y1_ub,"y0_lb"=y0_lb,"y0_ub"=y0_ub,
                 "ate_lb"=ate_lb,"ate_ub"=ate_ub,"ate_rps"=ate_rps)  

  class(outputs) = 'ATbounds'

  outputs
}

