#' @title Bounding the average treatment effect on the treated
#'
#' @description Bounds the average treatment effect on the treated under the unconfoundedness assumption without the overlap condition
#'
#' @param y n-dimensional vector of binary outcomes
#' @param t n-dimensional vector of binary treatments
#' @param x n by p matrix of covariates
#' @param rps n-dimensional vector of reference propensity scores
#' @param q polynomial order (default: q = 10)
#' @param discrete TRUE if x inclues only discrete covariates and FALSE if not (default: FALSE)
#' @return An S3 object of type "bter". The object has the following elements.
#' \item{lb}{the lower bound of ATT}
#' \item{ub}{the upper bound of ATT}
#'
#' @examples
#' # to be added
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
attbounds <- function(y, t, x, rps, q = 10L, discrete = FALSE){

  ### Nearest neighborhood estimation ###

  nn_data <- FNN::get.knnx(x, x, k=q)
  nn_i <- nn_data$nn.index
  nn_d <- nn_data$nn.dist

  nn_t <- {}
  for (k in 1:ncol(nn_i)){
    nn_t <- cbind(nn_t, t[nn_i[,k]])
  }

  if (discrete == TRUE){
    small_c <- 1e-8
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
    }
    if (q > 1){

      if ((q %% 2) == 1){ # if q is odd and q > 1
        v_x <- nx1 - nx1*(pxr^nx0)
        nx0c <- nx0 + (nx0 == 0L)
        rps_wt_nn <- -v_x/nx0c
      }

      if ((q %% 2) == 0){ # q is even
        v_x <- nx1 - nx*(pxr^nx0)
        nx0c <- nx0 + (nx0 == 0L)
        rps_wt_nn <- -v_x/nx0c
      }
    }

    att_lb <- t*(y-max(y)) + rps_wt_nn*(1-t)*(y-max(y))
    att_ub <- t*(y-min(y)) + rps_wt_nn*(1-t)*(y-min(y))
    att_lb <- sum(att_lb)/sum(t)
    att_ub <- sum(att_ub)/sum(t)

    outputs = list("lb"=att_lb,"ub"=att_ub)
    class(outputs) = 'ATbounds'

    outputs
}
