#' Generate a synthetic offspring
#'
#' Takes two haplotypes of equal length and create an offspring haplotype via recombination.
# 
# Arguments:
#' @param x1 One parental haplotype, a vector of length p with each entry 0 or 1.
#' @param x2 A second parental haplotype, a vector of length p with each entry 0 or 1.
#' @param d vector of genetic distances between sites. Should be the same length as x1. 
#' @param lambda recombination rate per unit of genetic distance
#
# Returns:
#' @return Offspring haplotype. A vector of length equal to x1 with each entry 0 or 1.
#' @export
generate_offspring = function(x1, x2, z = NULL, d = NULL, lambda = .012) {
  if(is.null(d)) {d = x1 * 0 + 1}
  if(is.null(z)) {
    z = generate_ancestry(d, lambda)
  }
  x = rep(0, length(x1))
  x[z == 0] = x1[z == 0]
  x[z == 1] = x2[z == 1]
  
  return(x)
}

#' Generate ancestry information
#'
#' Simulates recombination and returns vector of 0-1 indicating ancestry. 
#' Intended for intenal use with the \code{generate_offspring} function.
# 
# Arguments:
#' @param d Vector of genetic distances between sites. Should be the length of the haplotype.
#' @param lambda Recombination rate per unit of genetic distance.
#'
#' Return value:
#' @return  A vector of length equal to x1 with each entry 0 or 1, 
#'    corresponding to the ancestry at each location.
generate_ancestry = function(d, lambda = .012) {
  p = length(d)
  z = rep(0, p)
  
  u = runif(p) 
  recomb =  (u > (1 + exp(- 2 * d * lambda)) / 2) #recomb points
  z = cumsum(recomb) %% 2
  if(runif(1) < .5) {
    z = 1 - z
  }

  return(z)
}

#' Find a pair of haplotypes.
#'
#' Returns a pair of haplotypes. A first haplotype is uniformly selected, and then 
#'   a second is sampled with higher probability assigned to similar haplotypes.
#
# Arguments:
#' @param scores: A vector of the lenght of the parent population. Gives the score for 
#'    for each parent, such as the inner product with the 1st PC.
#
# Return values:
#' @return a vector of length two giving the selected parents.
find_mate_pair = function(scores, beta = 1) {
  i = sample(length(scores), 1)
  
  temp = scores[i] * scores
  probs = 1 / (1 + exp(-beta * temp))
  probs[i] = 0
  probs = probs / sum(probs)
  j = sample(length(scores), size = 1, prob = probs)
  
  return(c(i, j))
}

#' Simulate an offspring population
#' 
#' Samples a population of offspring haplotypes with optional preferred mating.
#
# Arguments:
#' @param parent_pop n x p matrix of parental haplotypes
#' @param d genetic distance
#' @param n_offspring number of offspring to sample.
#' @param lambda recombination rate (recombinations per unit genetic distance)
#' @param scores n-vector of scores for preferred mating (e.g. inner product with first PC)
#' @param beta intensity of preferred mating.
#
# Returns:
#' @return ancestry: n_offspring x 2 matrix giving parents for each sampled offspring
#' @return children: n_offspring x p matrix giving the offspring haplotypes.
create_offspring_population = function(parent_pop, d,  n_offspring = 100, 
                                       lambda = .012, scores = NULL, beta = 0) { 
  
  if(is.null(scores)) {
    scores = rep(1, nrow(parent_pop))
  }
  
  ancestry = replicate(n_offspring, find_mate_pair(scores, beta))
  ancestry = t(ancestry)
  
  children = apply(ancestry, 1, function(x){generate_offspring(parent_pop[x[1], ], parent_pop[x[2], ],
                                                                 d = d, lambda = lambda)})
  children = t(children)
  
  return(list(ancestry = ancestry, children = children))
}