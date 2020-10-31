#' @title Bounding the average treatment effect on the treated
#'
#' @description Bounds the average treatment effect on the treated under the unconfoundedness assumption without the overlap condition
#'
#' @param y n-dimensional vector of binary outcomes
#' @param t n-dimensional vector of binary treatments
#' @param x n by p matrix of covariates
#' @param rps n-dimensional vector of reference propensity scores
#' @param q polynomial order (default: q = 2, which uses the nearest neighbor excluding own observations)
#' @param permute_max maximum number of permutations to shuffle the data (default: 0)
#' @param discrete TRUE if x includes only discrete covariates and FALSE if not (default: FALSE)
#' @param studentize TRUE if x is studentized elementwise and FALSE if not (default: TRUE)
#' @param small_c a small positive constant to determine the two covariate vectors are identical (default: 1e-8)
#' This constatn is only used when the option 'discrete' is TRUE. 
#' 
#' @return An S3 object of type "ATbounds". The object has the following elements.
#' \item{lb}{the lower bound of ATT, i.e. E[Y(1) - Y(0) | T = 1]}
#' \item{ub}{the upper bound of ATT, i.e. E[Y(1) - Y(0) | T = 1]}
#' \item{att_rps}{the point estimate of ATT using the reference propensity score}
#' 
#' @examples
#'   Y <- RHC[,"survival"]
#'   D <- RHC[,"RHC"]
#'   X <- RHC[,-c(1,2)]
#'   rps <- rep(mean(D),length(D))  
#'   results_att <- attbounds(Y, D, X, rps, q = 3)
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
attbounds <- function(y, t, x, rps, q = 2L, permute_max = 0, discrete = FALSE, studentize = TRUE, small_c = 1e-8){
    
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
    
    # Reorder the sample to break ties in random permutation to shuffle the data 
    
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

  ### ATT estimation using reference propensity scores  ###    
        
  rps_wt <- rps/(1-rps)      
  att_rps <- sum(t*y-rps_wt*(1-t)*y)/sum(t)    
      
  ### Nearest neighbor estimation ###

  nn_data <- FNN::get.knnx(x, x, k=q)
  nn_i <- nn_data$nn.index
  nn_d <- nn_data$nn.dist

  nn_t <- {}
  for (k in 1:ncol(nn_i)){
    nn_t <- cbind(nn_t, t[nn_i[,k]])
  }

  if (discrete == TRUE){
    nx <- rowSums(nn_d < small_c)
    nx1 <- rowSums(nn_t*(nn_d < small_c))
  } else if (discrete == FALSE){
    nx <- q
    nx1 <- rowSums(nn_t)
  }
  else {
    stop("'discrete' must be either TRUE or FALSE.")
  }


  nx0 <- nx - nx1
  pxr <- rps/(rps-1)

  ### Computing weights and obtain bound estimates ###

    if (q == 1){
      rps_wt_nn <- 0
    } else if (q > 1){

      if ((q %% 2) == 1){ # if q is odd and q > 1
        v_x <- nx1 - nx1*(pxr^nx0)
      } else if ((q %% 2) == 0){ # q is even
        v_x <- nx1 - nx*(pxr^nx0)
      } else{
        stop("'q' must be a positive integer.")      
      }  
      
      nx0c <- nx0 + (nx0 == 0L)
      rps_wt_nn <- -v_x/nx0c
    } else{
      stop("'q' must be a positive integer.")      
    }

    att_lb <- t*(y-max(y)) + rps_wt_nn*(1-t)*(y-max(y))
    att_ub <- t*(y-min(y)) + rps_wt_nn*(1-t)*(y-min(y))
    att_lb <- sum(att_lb)/sum(t)
    att_ub <- sum(att_ub)/sum(t)

    result <- matrix(NA, nrow = 1, ncol = 2) 
    result <- c(att_lb, att_ub)
    results <- rbind(results,result)  
    
    }
    
    est <- apply(results,2,mean)
    
    if (FALSE){
    att_lb <- est[1]*(est[1] <= est[2]) + att_rps*(est[1] > est[2])
    att_ub <- est[2]*(est[1] <= est[2]) + att_rps*(est[1] > est[2])
    }
    
    att_lb <- est[1]
    att_ub <- est[2]
    
    outputs = list("lb"=att_lb,"ub"=att_ub,"att_rps"=att_rps)
    class(outputs) = 'ATbounds'

    outputs
}
