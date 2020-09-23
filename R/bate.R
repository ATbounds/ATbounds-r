#' @title Bounding the average treatment effect
#'
#' @description Bounds the average treatment effect under the unconfoundedness assumption without the overlap condition
#'
#' @param y n-dimensional vector of binary outcomes
#' @param t n-dimensional vector of binary treatments
#' @param x n by p matrix of covariates
#' @param rps n-dimensional vector of reference propensity scores
#' @param q polynomial order (default: q = 10)
#' @param discrete TRUE if x inclues only discrete covariates and FALSE if not (default: FALSE)
#' @return An S3 object of type "bter". The object has the following elements.
#' \item{y1_lb}{the lower bound of E[Y(1)]}
#' \item{y1_ub}{the upper bound of E[Y(1)]}
#' \item{y0_lb}{the lower bound of E[Y(0)]}
#' \item{y0_ub}{the upper bound of E[Y(0)]}
#' \item{ate_lb}{the lower bound of E[Y(1) - Y(0)]}
#' \item{att_ub}{the upper bound of E[Y(1) - Y(0)]}
#' 
#' @examples
#' # to be added
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
bate <- function(y, t, x, rps, q = 10L, discrete = FALSE){

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

  px1 <- ((rps-1)/rps)^nx1
  px0a <- ((rps-1)/rps)^nx0
  px0b <- (rps/(rps-1))^nx0

  ### Computing weights and obtain bound estimates ###

    if (q == 1){
      y1_wt <- 1
      y0_wt <- 1
    }
    if (q > 1){

      if ((q %% 2) == 1){ # if q is odd and q > 1
        v_x1 <- 1 - (nx0/nx)*px1
        v_x0 <- 1 - (nx1/nx)*px0b
        nx1c <- nx1 + (nx1 == 0L)
        nx0c <- nx0 + (nx0 == 0L)
        y1_wt <- nx*v_x1/nx1c
        y0_wt <- nx*v_x0/nx0c
      }

      if ((q %% 2) == 0){ # q is even
        v_x1 <- 1 - px1
        v_x0 <- 1 - px0a
        nx1c <- nx1 + (nx1 == 0L)
        nx0c <- nx0 + (nx0 == 0L)
        y1_wt <- nx*v_x1/nx1c
        y0_wt <- nx*v_x0/nx0c
      }
    }

    y1_lb <- min(y) + y1_wt*(treat == 1)*(y-min(y))
    y1_ub <- max(y) + y1_wt*(treat == 1)*(y-max(y))
    y1_lb <- mean(y1_lb)
    y1_ub <- mean(y1_ub)
    
    y0_lb <- min(y) + y0_wt*(treat == 0)*(y-min(y))
    y0_ub <- max(y) + y0_wt*(treat == 0)*(y-max(y))
    y0_lb <- mean(y0_lb)
    y0_ub <- mean(y0_ub)
    
   ate_lb <-  y1_lb - y0_ub
   ate_ub <-  y1_ub - y0_lb
   

    outputs = list("y1_lb"=y1_lb,"y1_ub"=y1_ub,
                   "y0_lb"=y0_lb,"y0_ub"=y0_ub,
                   "ate_lb"=ate_lb,"ate_ub"=ate_ub)
    class(outputs) = 'batt'

    outputs
}

# this code is really for continuous covariates
# fir discrete covariates, we find unique elements and 
# then estimate the empirical prob. mass function.
