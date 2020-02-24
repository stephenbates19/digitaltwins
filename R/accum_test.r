#' Hinge-exponential function 
#'
#' Intended for internal use with the \code{accum_test} function.
# Arguments:
#' @param x A vector of p-values.
#' @param c A parameter of the accumulation test.
# Returns:
#' @return A vector with the same length as x of 
#'    the element-wise application of 
#' 		the hinge-exponential function.
hinge_exp <- function(x, c = 2) {
  out <- c * log(1 / (c * (1 - x)))
  out[out <= 0 ] <- 0
  return(out)
}

#' Accumulation test with the hinge-exponential function
#'
#' Takes an ordered list of p-values and returns a set of selections
#' controlling the false discovery rate.
#' This is an implementation of the method of Li and Barber, JASA 2017. 
# Arguments:
#' @param x A vector of p-values to be tested IN ORDER.
#' @param	alpha The desired FDR threshold.
#' @param c Parameter of the hinge-exponential function used in the test.
#' @param strict If TRUE, the procedure is more conservative but
#'  controls the false discovery rate. If FALSE, the procedure is more liberal
#'  but controls a modified version of the false discovery rate.
# Returns:
#' @return A vector of rejections, has the form 1:k for some value of k.
#' @export
accum_test <- function(x, alpha = .2, c = 2, strict = FALSE) {
  if(strict) {
    accum_val <- (c+cumsum(hinge_exp(x, c))) / seq(2, 1+length(x))
  } else {
    accum_val <- cumsum(hinge_exp(x, c)) / seq(1, length(x))
  }
  below <- which(accum_val <= alpha)
  if(length(below)==0) {
    return(c())
  } else {
    return(1:max(below))
  }
}
