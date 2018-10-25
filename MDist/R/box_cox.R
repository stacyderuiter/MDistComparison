#' box-cox and inverse box-cox transform
#'
#' Box-Cox and inverse Box-Cox power transformation (often used to normalize).
#' @param x data
#' @param lambda transform parameter
#' @param inverse logical.  Inverse transform? Defaults is FALSE (use inverse=TRUE to reverse the transform, i.e. back to natural units)
#' @return a matrix of transformed data that is as close to multi-variate normal as possible with this transform; different variables are in different columns
#' @export

box_cox <- function(x, lambda, inverse=FALSE) {
  #x should be a matrix of dive phase durations
  #lmat should be a matrix r of transform params
  if (!inverse){
    (x^lambda-1)/lambda
  }else{
    (x*lambda + 1)^(1/lambda)
  }
}
