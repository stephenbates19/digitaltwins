# sample_parity_poisson
#' Sample poisson conditional on its parity
#'
#' Sample a poisson conditional on it being even or odd
#'   This function is only computationally efficient when lambda is small.
#
# Arguments:
#' @param lambda Poisson parameter.
#' @param parity Either "even" or "odd".
#
# Returns:
#' @return Scalar, a sample from the poisson conditional on the specified parity
sample_parity_poisson <- function(lambda, parity = "even") {
  ##choose starting value and normalizing constant
  if(parity == "even") {
    i <- 0
    alpha <- .5 * (1 + exp(-2 * lambda)) ##normalizing constant
  } else if(parity == "odd") {
    i <- 1
    alpha <- 1 - .5 * (1 + exp(-2 * lambda)) ##normalizing constant
  } else {
    warning("invalid parity specified")
    return(0)
  }

  u <- runif(1)
  cur_mass <- exp(-lambda) * lambda^(i) / factorial(i) / alpha
  ##loop until CDF is greater than the random uniform variable for the first time
  while(cur_mass <= u & i < 100) {
    i <- i + 2
    cur_mass  <- cur_mass + exp(-lambda) * lambda^(i) / factorial(i) / alpha
  }
  ##if(i == 100) {print(u)}

  return(i)
}

# impute_variant_fast
#' Imputes a set of variants
#'
#' Imputes the variants in a specified group, conditional on other variants,
#'  based on the recombination model. This is intented for use a subroutine of
#'   the conditional randomization test.
#
# Arguments:
#' @param x 0-1 vector of length p. The observed haplotypes of the offspring.
#' @param p1 0-1 vector of length p. One of the parental haplotypes 
#'    giving rise to \code{x}.
#' @param p2 0-1 vector of length p. The other parental haplotypes 
#'    giving rise to \code{x}.
#' @param group a subset of indices of 1:p. This subset of variants will be imputed.
#' @param d numeric vector of length p. The genetic distance between variants
#' @param lambda the recombination rate, in units of the vector d
#' @param epsilon the pointwise mutation rate
#' @param window_size when imputing the variants, look only at a window of this size
#'     surrounding the group g. A smaller window size leads to faster computation
#'     and avoids numerical issues.
#
# returns:
#' @return A vector of the same length as \code{group} containing the imputed variants.
#' @export
impute_variant_fast <- function(x, p1, p2, group, d, fb_result = NULL,
                                lambda = .012, epsilon = .001, window_size = 200) {
  start <- group[1]
  end <- group[length(group)]
  p <- length(x)

  ## check case where group is entire chromosome,
  if(length(x) == length(group)) {
    return(generate_offspring(p1, p2, d = d, lambda = lambda))
  }

  ##to avoid numerical innacuracies, run FB only on a subset of chromosome near the group
  start_window <- max(1, start - window_size)
  end_window <- min(p, end + window_size)
  start_j <- start - start_window ##index of start of group, 0-indexed
  ##check if fb_result has been pre-computed and saved
  if(is.null(fb_result)) {
    window <- start_window:end_window
    fb_result <- digitaltwins:::forward_backward_(x[window], p1[window], p2[window], d[window], length(x[window]),
                                                lambda, epsilon, start_j, start_j + length(group) - 1)
  }

  ##compute probability of and odd number of recombinations
  p1_start <- fb_result[1]
  z1 <- (runif(1) < p1_start) ##starting state
  if(z1 == TRUE) {
    p_odd <- fb_result[2]
  } else {
    p_odd <- fb_result[3]
  }

  ##sample recombination events
  d_total <- sum(d[group])
  if(runif(1) < p_odd) {
    n_recomb <- sample_parity_poisson(d_total * lambda, parity = "odd")
  } else {
    n_recomb <- sample_parity_poisson(d_total * lambda, parity = "even")
  }
  if(n_recomb == 0) {
    if(z1 == TRUE) {
      out <- p1[group]
    } else {
      out <- p2[group]
    }
    corruptions <- rbinom(n = length(group), size = 1, prob = epsilon)
    out[corruptions == 1] <- 1 - out[corruptions == 1]
    return(out)
  } else if(n_recomb >= length(group)) {
    recomb_points <- group
  } else {
    prob <- d[group] / sum(d[group])
    recomb_points <- sample(length(group), size = n_recomb, replace = FALSE, prob = prob)
  }

  ##construct the CRT replicate
  z <- rep(0, length(group) + 1) ##padding on left end
  z[1] <- as.integer(z1==FALSE)
  z[recomb_points+1] <- 1
  z <- (cumsum(z) %% 2)
  z_no_pad <- z[2:(length(group) + 1)]
  x_prime <- p1[group]
  x_prime[z_no_pad == 1] <- p2[group[z_no_pad == 1]]

  ##add transcription errors
  corruptions <- rbinom(n = length(group), size = 1, prob = epsilon)
  x_prime[corruptions == 1] <- 1 - x_prime[corruptions == 1]
  return(x_prime)
}

# get_fb_results
#' Forward-backward algorithm for parent-offspring trio of haplotypes
#'
#' Returns the result of the forward backward algorithm for a sepcified group.
#'   Intended for use as a pre-computation step for the \code{linear_crt} function.
#
# Arguments:
#' @param H an (n x p) matrix of the subject haplotypes
#' @param H_parent and (n_2 x p) matrix of the parental haplotypes
#' @param anc a (n x 2) table of ancestries.
#'     anc[i, 1] and and[i, 2] give the rows of H_parent corresponding to the
#'     parents of subject i
#' @param   group a list of of indices of the group to test. Should be a continuous region,
#'     such as 10:20.
#' @param   d a vector of length p giving the genetic distances between sites
#' @param   lambda the recombination rate, in units of d
#' @param   epsilon the pointwise mutation rate
#' @param   window_size the number of sites to on each size of the group to
#'    use when computing the forward-backward probabilities. Smaller groups improve speed
#'    and avoid numerical issues.
#
# Returns:
#' @return a (n x 3) matrix giving the results of the forward-backward algorithm.
#' See the documentation of the \code{forward_backward_} function for output information.
get_fb_results <- function(H, H_parent, anc, group, d, lambda = .012, epsilon = .001, window_size = 200) {
  start <- group[1]
  end <- group[length(group)]
  p <- ncol(H)

  ## Define window
  start_window <- max(1, start - window_size)
  end_window <- min(p, end + window_size)
  start_j <- start - start_window ##index of start of group, 0-indexed

  FB.WORKER <- function(i) {
    digitaltwins:::forward_backward_(H[i, start_window:end_window],
                                   H_parent[anc[i, 1], start_window:end_window],
                                   H_parent[anc[i, 2], start_window:end_window],
                                   d[start_window:end_window],
                                   length(H[i, start_window:end_window]),
                                   lambda, epsilon, start_j, start_j + length(group) - 1)
  }
  output <- t(vapply(1:nrow(H), FB.WORKER, double(3)))
  return(output)
}
