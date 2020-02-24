#' Computes TDT p-value
#'
#' Arguments:
#' @param y 0-1 vector, outcome status
#' @param child 0-1 vector, haplotype of child
#' @param p1 0-1 vector, haplotype of one parent
#' @param p2 0-1 vector, haplotype of second parent
#
# returns:
#' @return p-value
#' @export
tdt = function(y, child, p1, p2) {
  a = sum((y == 1) & (child == 1) & (p1 != p2))
  b = sum((y == 1) & (child == 0) & (p1 != p2))
  chisq_stat = (a - b)^2 / (a + b)
  
  1 - pchisq(chisq_stat, df = 1)
}



