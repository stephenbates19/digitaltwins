#include <Rcpp.h>
using namespace Rcpp;

//' forward_backward_
//'
//' Carry out the forward-backward algorithm for the HMM implied by knowing the two
//' parental haplotypes.
//'
//' @param x  0-1 vector of length p giving the offspring haplotypes
//' @param p1 0-1 vector of length p giving one parental haplotype
//' @param p2 0-1 vector of length p giving the other parental haplotype
//' @param d numberic vector of length p - 1  giving the distance between each site
//' @param p length of the vector x
//' @param lambda the recombination rate, a multiplier for the distance d
//' @param epsilon the pointwise mutation rate
//' @param start indexes of the subset of x to perform the computation (1 indexed)
//' @param end indexes of the subset of x to perform the computation (1 indexed)
//'
//' @return p1_start: the probability that the latent state immediately before the specified window is 1,
//'        i.e. the offspring is being copied from p1.
//' @return p1_start_odd: probability of an odd number of recombinations given that 
//'       the latent state immediately before the specified window is 1.
//' @return  p2_start_odd: probability of an odd number of recombinations given that 
//'       the latent state immediately before the specified window is 2.
// [[Rcpp::export]]
NumericVector forward_backward_(NumericVector x, NumericVector p1, 
                               NumericVector p2, NumericVector d, int p, double lambda, double epsilon,
                               int start, int end) {

  NumericMatrix a(p, 2); //forward probs
  NumericMatrix b(p, 2); //backward probs
  NumericVector out(3); //store return values
  
  //forward pass
  double p_tran;
  a(0, 0) = (x[0] == p1[0]) * (1 - epsilon) + (x[0] != p1[0]) * epsilon;
  a(0, 1) = (x[0] == p2[0]) * (1 - epsilon) + (x[0] != p2[0]) * epsilon;
  for(int j = 1; j < p; j++) {
    p_tran = 1 - .5 * (1 + exp(-2 * d(j) * lambda));
    a(j, 0) = (a(j-1, 0) * (1 - p_tran)  + a(j-1, 1) * p_tran)* ((x[j] == p1[j]) * (1 - epsilon) + (x[j] != p1[j]) * epsilon);
    a(j, 1) = (a(j-1, 1) * (1 - p_tran)  + a(j-1, 0) * p_tran)* ((x[j] == p2[j]) * (1 - epsilon) + (x[j] != p2[j]) * epsilon);
  }
  
  //backward pass
  b(p-1, 0) = 1;
  b(p-1, 1) = 1;
  for(int j = p-2; j >= 0; j--) {
    p_tran = 1 - .5 * (1 + exp(-2 * d(j) * lambda));
    b(j, 0) = b(j + 1, 0) * (1 - p_tran) * ((x[j+1] == p1[j+1]) * (1 - epsilon) + (x[j+1] != p1[j+1]) * epsilon) + (
      b(j + 1, 1) * p_tran * ((x[j+1] == p2[j+1]) * (1 - epsilon) + (x[j+1] != p2[j+1]) * epsilon));
    b(j, 1) = b(j + 1, 0) *  p_tran * ((x[j+1] == p1[j+1]) * (1 - epsilon) + (x[j+1] != p1[j+1]) * epsilon) + (
      b(j + 1, 1) * (1 - p_tran) * ((x[j+1] == p2[j+1]) * (1 - epsilon) + (x[j+1] != p2[j+1]) * epsilon));
  }
    
  //extract results
  // compute probability of an odd number of transitions in the window
  double d_total = 0;
  for(int j = start; j <= end; j++) {
    d_total += d[j];
  }
  double p_even_uncond = .5 * (1 + exp(-2 * lambda * d_total));
  double p1_start_odd; // probability of odd # recombinations, if starting in state 1
  double p2_start_odd; // probability of odd # of recombinations, if starting in state 2
  if(start == 0 || end == p - 1) {
    p1_start_odd = (1 - p_even_uncond);
    p2_start_odd = (1 - p_even_uncond);
  } else {
    p1_start_odd = (1 - p_even_uncond) * (a(start - 1, 0) * b(end, 1));
    double p1_start_even = p_even_uncond * (a(start - 1, 0) * b(end, 0));
    p1_start_odd = p1_start_odd / (p1_start_odd + p1_start_even);
    p2_start_odd = (1 - p_even_uncond) * (a(start - 1, 1) * b(end, 0));
    double p2_start_even = p_even_uncond * (a(start - 1, 1) * b(end, 1));
    p2_start_odd = p2_start_odd / (p2_start_odd + p2_start_even);
  }
  //compute probability dist of latent variable adjacent to the start of the group
  double p1_start; 
  if(start == 0) {
    p1_start = .5;
  } else if (end == p - 1) {
    p1_start = a(start - 1, 0) / (a(start - 1, 0) + a(start - 1, 1));
  } else {
    p1_start = a(start-1, 0) * b(end, 0) * p_even_uncond + a(start-1, 0) * b(end, 1) * (1 - p_even_uncond); 
    double p2_start = a(start-1, 1) * b(end, 1) * p_even_uncond + a(start, 1) * b(end, 0) * (1 - p_even_uncond); 
    p1_start = p1_start / (p1_start + p2_start);
  }

    
  out[0] = p1_start;
  out[1] = p1_start_odd;
  out[2] = p2_start_odd;
  // out[3] = a(start - 1, 0);
  // out[4] = a(start - 1, 1);
  // out[5] = b(end, 0);
  // out[6] = b(end, 1);
  
  return(out);
}


//'sample_ancestry_
//'
//' Sample the latent state of ancenstry indicators for a haplotype given
//'   its two ancestral haplotypes.
//'
//' @param x 0-1 vector of length p giving the offspring haplotypes
//' @param p1 0-1 vector of length p giving one parental haplotype
//' @param p2 0-1 vector of length p giving the other parental haplotype
//' @param d numberic vector of length p - 1  giving the distance between each site
//' @param p length of the vector x
//' @param lambda the recombination rate, a multiplier for the distance d
//' @param epsilon the pointwise mutation rate
//'
//' @return a vector of of the ancestral states
// [[Rcpp::export]]
NumericVector sample_ancestry_(NumericVector x, NumericVector p1, 
                               NumericVector p2, NumericVector d, int p, double lambda, double epsilon) {

  NumericMatrix b(p, 2); //backward probs
  NumericVector out(p); //ancestral states
  double p_tran;
  double p_tot;

  //backward pass
  b(p-1, 0) = 1;
  b(p-1, 1) = 1;
  for(int j = p-2; j >= 0; j--) {
    p_tran = 1 - .5 * (1 + exp(-2 * d(j) * lambda));
    b(j, 0) = b(j + 1, 0) * (1 - p_tran) * ((x[j+1] == p1[j+1]) * (1 - epsilon) + (x[j+1] != p1[j+1]) * epsilon) + (
      b(j + 1, 1) * p_tran * ((x[j+1] == p2[j+1]) * (1 - epsilon) + (x[j+1] != p2[j+1]) * epsilon));
    b(j, 1) = b(j + 1, 0) *  p_tran * ((x[j+1] == p1[j+1]) * (1 - epsilon) + (x[j+1] != p1[j+1]) * epsilon) + (
      b(j + 1, 1) * (1 - p_tran) * ((x[j+1] == p2[j+1]) * (1 - epsilon) + (x[j+1] != p2[j+1]) * epsilon));
  }

  // random number generation
  RNGScope scope;
  NumericVector u = runif(p);

  // sample the Markov chain
  for(int j = 0; j < p; j++) {
    // terms from node j and j+1:p
    double p0 = (x[j] == p1[j]) * (1 - epsilon) + (x[j] != p1[j]) * epsilon;
    p0 = p0 * b(j, 0);
    double p1 = (x[j] == p2[j]) * (1 - epsilon) + (x[j] != p2[j]) * epsilon;
    p1 = p1 * b(j, 1);
    
    // terms from previous state
    p_tran = 1 - .5 * (1 + exp(-2 * d(j) * lambda));
    if(j != 0) {
      if(out(j-1) == 0) {
        p0 = p0 * (1 - p_tran);
        p1 = p1 * p_tran;
      } else {
        p0 = p0 * p_tran;
        p1 = p1 * (1 - p_tran);
      }
    }

    // compute combined probability and sample
    p_tot = p0 / (p0 + p1);
    if(u(j) < p_tot) {
      out(j) = 0;
    } else {
      out(j) = 1;
    }
  }

  return(out);
}
