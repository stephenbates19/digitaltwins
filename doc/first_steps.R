## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(causalsnps)

## ------------------------------------------------------------------------
#simulation parmaters
n_ext <- 1000 #number of external observations
n_trio <- 250 #number of trios
p <- 500 #number of observed genetic variants
chrome_width <- 50 #length of chromosome
d <- rep(p / chrome_width, p) #genetic distance between sites

k <- 5 #number of causal variants

## ------------------------------------------------------------------------
 data("haps_matrix")

## ------------------------------------------------------------------------
external_haps <- haps_matrix[1:(2*n_ext), ]
dim(external_haps)

## ------------------------------------------------------------------------
parent_haps <- haps_matrix[(2*n_ext + 1):(2*n_ext + 4 * n_trio), ]
dim(parent_haps)

## ------------------------------------------------------------------------
anc <- matrix(1:(4*n_trio), nrow = 2*n_trio, byrow = TRUE) #index of ancestors of haplotypes
print(dim(anc))
head(anc)

## ------------------------------------------------------------------------
set.seed(300)
# use the "generate_offspring" function for each row of the ancestry table
offspring_haps <- mapply(
  function(i, j) {generate_offspring(parent_haps[i, ], 
                                     parent_haps[i, ], 
                                     d = d)},
  anc[, 1], anc[, 2])
dim(offspring_haps)

## ------------------------------------------------------------------------
#create the regression coefficients
beta <- rep(0, p)
causal_variants <- sample(1:p, k)
beta[causal_variants] <- 1

#sample the response variable from a sparse linear model
Y_ext <- haps_to_gen(external_haps) %*% beta + rnorm(n_ext)
Y_offspring <- haps_to_gen(offspring_haps) %*% beta + rnorm(n_trio)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(glmnet)

## ------------------------------------------------------------------------
lasso_fit <- cv.glmnet(haps_to_gen(external_haps), Y_ext)
beta_hat <- coef(lasso_fit)
length(beta_hat) #Fitted linear predictor. First entry is an intercept.

## ------------------------------------------------------------------------
test_region <- 1:100
sum(causal_variants %in% test_region) #number of causal variants in the test region

## ------------------------------------------------------------------------
p_value <- linear_crt(offspring_haps, parent_haps, anc, Y_offspring, 
                      matrix(beta_hat, ncol = 1), 
                      group = test_region, d = d, family = "gaussian")
#null p-value
p_value 

## ------------------------------------------------------------------------
test_region <- 400:500
sum(causal_variants %in% test_region) #number of causal variants in the test region

## ------------------------------------------------------------------------
p_value <- linear_crt(offspring_haps, parent_haps, anc, Y_offspring, 
                      matrix(beta_hat, ncol = 1), 
                      group = test_region, d = d, family = "gaussian")

#non-null p-value
p_value 

