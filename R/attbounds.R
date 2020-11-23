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
attbounds <- function(Y, D, X, rps, Q = 2L, n_permute = 0, studentize = TRUE){
    
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
        att_wt <- 0
    
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
      pxr <- rps/(rps-1)

      if ((Q %% 2) == 1){ # if Q is odd and Q > 1
        v_x <- nx1 - nx1*(pxr^nx0)
      } else if ((Q %% 2) == 0){ # Q is even
        v_x <- nx1 - nx*(pxr^nx0)
      } else{
        stop("'Q' must be a positive integer.")      
      }  
      
      nx0c <- nx0 + (nx0 == 0L)
      att_wt <- v_x/nx0c
      
    } else{
      stop("'Q' must be a positive integer.")      
    }

    att_lb <- D*(Y-max(Y)) - att_wt*(1-D)*(Y-max(Y))
    att_ub <- D*(Y-min(Y)) - att_wt*(1-D)*(Y-min(Y))
    att_lb <- sum(att_lb)/sum(D)
    att_ub <- sum(att_ub)/sum(D)

    result <- matrix(NA, nrow = 1, ncol = 2) 
    result <- c(att_lb, att_ub)
    results <- rbind(results,result)  
    
    }
    
    est <- apply(results,2,mean)
    
    att_lb <- est[1]
    att_ub <- est[2]
    
    outputs = list("lb"=att_lb,"ub"=att_ub,"att_rps"=att_rps)
    class(outputs) = 'ATbounds'

    outputs
}
