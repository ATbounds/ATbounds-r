#' @title Summary method for ATbounds objects
#'
#' @description Produce a summary for an ATbounds object.

#' @param object ATbounds object
#' @param ... Additional arguments for summary generic
#' 
#' @return A summary is produced with bounds estimates and confidence intervals.
#' In addition, it has the following elements.
#' \item{Lower_Bound}{lower bound estimate and lower end point of the confidence interval}
#' \item{Upper_Bound}{upper bound estimate and upper end point of the confidence interval}
#' 
#' @examples
#'   Y <- RHC[,"survival"]
#'   D <- RHC[,"RHC"]
#'   X <- RHC[,c("age","edu")]
#'   rps <- rep(mean(D),length(D))
#'   results_ate <- atebounds(Y, D, X, rps, Q = 3)
#'   summary(results_ate)
#'
#' @references Sokbae Lee and Martin Weidner. Bounding Treatment Effects by Pooling Limited Information across Observations.
#'
#' @export
summary.ATbounds <- function(object,...){

  heading=c(paste("ATbounds: ",object$type,sep=""),
            paste("Call:",format(object$call)),
            paste("Confidence Level:", format(object$cov_prob)))
  
  est_lb = object$est_lb
  est_ub = object$est_ub
  ci_lb = object$ci_lb
  ci_ub = object$ci_ub
  
  sumob = data.frame(Lower_Bound = c(est_lb,ci_lb), Upper_Bound = c(est_ub,ci_ub), 
                     row.names = c("Bound Estimate", "Confidence Interval"))
  
  structure(sumob,heading=heading,class=c("anova","data.frame"))
  
}
