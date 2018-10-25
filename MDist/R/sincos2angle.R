#' convert sine and cosine to angle
#'
#' convert sine and cosine to an angle
#' @param sin sine of the angle
#' @param cos cosine of the angle
#' @return the angle
#' @export

sincos2angle <- function(sin, cos) {
  i <- complex(real = 0, imaginary = 1)
  ang <- Arg(cos + i * sin)
  return(ang)
  # three cheers for Euler!
}