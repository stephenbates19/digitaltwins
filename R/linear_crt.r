#' Run CV glmnet
#'
#' Returns the fitted coefficients from a lasso fit after tuning the
#'   regularization parameter by cross-validation. Intended for internal
#'   use with the \code{linear_crt} function.
#
# Arguments:
#' @param X and n by p matrix
#' @param Y a response vector of length n
#' @param train_idx a set of indices to train on
#' @param family type of regression, either "gaussian" or "binomial"
#' @param s coefficient vector to  return. either "lambda.min" or "lambda.1se"
#' @param parallel set to TRUE to run CV in parallel. Must initialize DoMC or DoParallel
#'     prior to execution.
#
# Returns:
#' @return Fitted coefficients, a vector of length p + 1.
#'      First coordinate is the intercept.
get_beta_glmnet <- function(X, Y, train_idx, family = "gaussian", s = "lambda.min", parallel = FALSE) {
  cv_fit <- cv.glmnet(X[train_idx, ], Y[train_idx], family = family, nfolds = 5, parallel = parallel)
  return(coef(cv_fit, s = s))
}

#' Logistic-log likelihood
#'
#' Returns the log-likelihood of logistic regression fit(s). Intended for use
#'   as a feature-importance measure for the CRT.
#
# Arguments:
#' @param logit_inv Aector of length n OR a matrix of dimension (n x k) giving the logit inverse of the probability
#'     of being one for k different models.
#' @param Y Observed responses, a 0-1 vector of length n or a matrix of dimension (n x k).
#' @param genotypes Whether to interpret the data as genotypes, i.e. summing pairs of rows
#'    of \code{logit_inv} gives the total logit inverse for the observation. 
#'
# Returns:
#' @return A vector of length k, the log-likelihoods of the columns of Y
logistic_ll <- function(logit_inv, Y, genotypes = TRUE) {
  if(genotypes) {
    return(logistic_ll(haps_to_gen(logit_inv), Y, genotypes = FALSE))
  }
  output <- sapply(1:ncol(Y), function(i) {
    p_vec <- 1 / (1 + exp(-logit_inv[,i]))
    out <- sum((Y[, i] == 0) * log(1 - p_vec) + (Y[, i] == 1) * log(p_vec))
    return(out)
  })
  return(output)
}

#' Correlation across columns
#'
#' Returns correlation of a prediction direction(s) with an observed response(s).
#'   Intended for use as a feature-importance measure for the CRT.
#
# Arguments:
#' @param  preds Vector of length n OR a matrix of dimension (n x k) giving the logit inverse of the probability
#'     of being one for k different models.
#' @param  Y Observed responses, a 0-1 vector of length n or a matrix of dimension (n x k)
#' @param genotypes Whether to interpret the data as genotypes, i.e. summing pairs of rows
#'    of \code{logit_inv} gives the total logit inverse for the observation. 
#'
#' Returns:
#' @return A vector of length k, the correlations of the columns of Y
col_cors <- function(preds, Y, genotypes = TRUE) {
  if(genotypes) {
    return(col_cors(haps_to_gen(preds), Y, genotypes = FALSE))
  }
  output <- sapply(1:ncol(Y), function(i) cor(preds[,i], Y[,i]))
  return(output)
}

#' Convert haplotype matrix to genotype matrix
#'
#' Adds pairs of rows of a matrix to convert a haplotype matrix
#' to a genotype matrix. For example, rows 1 and 2 are added to create
#' row 1 of the output matrix. Rows 3 and 4 are added to create row 2
#' of the output matrix, and so on.
#'
#' This function can also handle general arrays. It assumes that the 
#' first dimension indexes the haplotypes.
#'
# Argmuents:
#' @param H An array with an even length in the first dimension.
#'
#' @return An array of the same dimension as H, with the first dimension half as long.
#' @export
haps_to_gen <- function(H) {
  n <- nrow(H) / 2
  return(H[(1:n)*2,,drop=F] + H[(1:n)*2-1,,drop=F])
}


#' Complete haplotype list
#'
#' Returns a vector of all haplotypes from a set of individuals, given
#'   a list of some of the haplotypes from those individuals.
#'
# Arguments:
#' @param ids Vector of ids (integers)
#
# returns:
#' @return Vector of ids, typically of size 2 * length(ids)
get_genotype_ids <- function(ids) {
  out <- c( ((ids+1) %/% 2) * 2, ((ids+1) %/% 2) * 2 - 1)
  return(sort(unique(out)))
}

#' Convert genotype IDs to haplotype IDs
#'
# Arguments:
#' @param ids vector of ids (integers)
#'
# Returns:
#' @return Vector of ids, typically of size 2 * length(ids).
gen_id_to_haps <- function(ids) {
  return(sort(c(2*ids, 2*ids - 1)))
}


#' Convert haplotype IDs to genptype IDs
#'
# Arguments:
#' @param  ids vector of ids (integers)
#
# Returns:
#' @return Vector of ids, typically of size length(ids) / 2
haps_id_to_gen <- function(ids) {
  return(sort(unique((ids - 1) %/% 2 +  1)))
}

#' Computed test statistics based on deviance change for GLM
#'
# Arguments:
#' @param  pred_base vector of baseline predictions
#' @param  x matrix of new predictors
#' @param  y vector of response observations
#
# Returns:
#' @return A test statistic
deviance_test_stats <- function(pred_base, x, y, family="gaussian") {
  model.1 <- glm(y~offset(pred_base), family=family)
  model.2 <- glm(y~pred_base+x, family=family)
  model.comparison <- anova(model.1,model.2,test="Chisq")
  if(length(model.comparison$Deviance)<0) return(0)
  deviance <- model.comparison$Deviance[2]
  if(is.na(deviance)) return(0)
  return(deviance)
}


#' DTT using a fixed linear statistic
#'
#' Returns p-values from the DTT using a linear feature importance measure.
#
# Arguments:
#' @param H an (2n x p) matrix of the subject haplotypes.
#'     It is assumed rows 1,2 belong to subject 1, rows 3,4 belong to subject 2, etc.
#' @param H_parent and (n_2 x p) matrix of the parental haplotypes
#' @param anc a (2n x 2) table of ancestries.
#'     anc[i, 1] and and[i, 2] give the rows of H_parent corresponding to the
#'     parents of row i of the haplotypes H
#' @param Y observed responses, a vector of length n or a matrix of dimension (n x k)
#' @param beta feature importance directions, fit on from an ind data set.
#'     A vector of length p+1 or a matrix of dimension ((p+1) x k).
#'     e.g. The output of get_beta_glment. The first coefficient is assumed to
#'     be an intercept.
#' @param d a vector of length p of genentic distances
#' @param adjust a vector of length n or a matrix of dimension (n x k),
#'     giving the contribution of the other chromosomes to the likelihood.
#'     I.e. adjust = X %*% beta_hat where X is all other chromosomes and beta is
#'     the fitted coefficients.
#' @param group a list of of indices of the group to test. Should be a continuous region,
#'     such as 10:20.
#' @param test_idx a set of indices (between 1 and 2n) used to compute the test statistics
#' @param family type of regression. Either "gaussian" or "binomial".
#'     If "guassian", correlation is used as a feature importance statistic. If "binomial",
#'    logistic log-likelihood is used as a feature importance statistic.
#' @param n_reps number of repitions of the CRT to carry out
#' @param genotypes Defaults to TRUE. If FALSE, all rows are assumed to be independent
#'     (i.e. no two haplotypes are from the same individual).
#' @param verbose if TRUE, prints various diagnostic messages to console
#' @param parallel requires doMC to be registered (default: FALSE)
#' @param recomb_thresh threshold of probability of recombination events to ignore. Lower
#'	   values will have potentially higher power at increased computational cost.
#' @param f_stat whether or not to use the f-statistic for continuous responses.
#'	   The f-statistic may be much slower for groups with many nonzero coefficients.
#
# Returns:
#' @return A vector of length k, a p-value for each column of Y.
#' @export
linear_crt <- function(H, H_parent, anc, Y, beta, group, d, adjust = NULL,
                       test_idx = 1:nrow(H), family = "gaussian", n_reps = 100,
                       genotypes = TRUE, verbose = FALSE, parallel = FALSE,
                       recomb_thresh = 0.001, f_stat = FALSE) {

  stopifnot( is.matrix(H) )
  stopifnot( is.matrix(H_parent) )
  stopifnot( is.matrix(anc) )

  ## Check dimensions of input
  stopifnot( ncol(H) == ncol(H_parent) )
  #stopifnot( (2*nrow(H)) == nrow(H_parent) )
  stopifnot( nrow(H) == nrow(anc) )
  stopifnot( ncol(anc) == 2 )
  stopifnot( length(group) <= ncol(H) )
  stopifnot( length(d) == ncol(H) )
  stopifnot( length(test_idx) <= nrow(H) )

  if(is.null(dim(beta))) {
    beta <- as.matrix(beta, ncol = 1)
  }
  if(is.null(dim(Y))) {
    Y <- as.matrix(Y, nrow=length(Y), ncol = 1)
  }
  stopifnot( ncol(Y) == ncol(beta) )
  stopifnot( (2*nrow(Y)) == nrow(H) )
  stopifnot( nrow(beta) == (ncol(H)+1) )

  if(is.null(adjust)) {
    adjust <- matrix(0,nrow(Y),ncol(beta))
  }
  if(is.null(dim(adjust))) {
    adjust <- as.matrix(adjust,length(adjust),1)
  }
  stopifnot( nrow(adjust) == nrow(Y) )
  stopifnot( ncol(adjust) == ncol(Y) )

  if(ncol(H) != length(group)) {
    ## Case 1: group is not entire chromosome
    ## Pre-compute forward-backward results
    cat("Pre-computing forward-backward results... ")
    fb_results <- get_fb_results(H[test_idx, ], H_parent, anc[test_idx, ], group, d,
                                 lambda = .012, epsilon = .001, window_size = 200)
    p_odd <- fb_results[, 1] * fb_results[, 2] + (1 - fb_results[, 1]) * fb_results[, 3]
    fb_full <- matrix(0, nrow = nrow(H), ncol = 3)
    fb_full[test_idx, ] <- fb_results
    cat("done.\n")

    ## Ignore observations with no variability (Matteo: this would decrease power)
    ## The total number of recombinations will be much larger if we ignore none
    recomb_idx <- which(p_odd >= recomb_thresh)
    impute_idx <- test_idx[recomb_idx]
    if(length(impute_idx) == 0) {
      cat("Warning, no recombination events detected.\n")
      return(rep(1, ncol(beta)))
    } else {
      cat(sprintf("Number of possible recombination events detected: %d.\n", length(impute_idx)))
    }
  } else {
    ## Case 2: group is entire chromosome
    impute_idx <- test_idx
  }

  ## Select ids to use in test statistic
  if(genotypes) {
    eval_idx <- get_genotype_ids(impute_idx)
    eval_y <- haps_id_to_gen(eval_idx)
  } else {
    eval_idx <- impute_idx
    eval_y <- eval_idx
  }

  cat("Computing stats... ")
  ## Pre-compute some test stat information
  if(genotypes) {
    intercept <- matrix(beta[1,,drop=F], nrow=length(eval_idx), ncol=ncol(beta), byrow=TRUE) / 2
    adjust_mat <- adjust[rep(1:nrow(adjust),each=2),,drop=F][eval_idx,,drop=F] / 2
  } else {
    intercept <- matrix(beta[1,,drop=F], nrow=length(eval_idx), ncol=ncol(beta), byrow=TRUE)
    adjust_mat <- adjust[eval_idx,,drop=F]
  }
  pred_base <- intercept + adjust_mat ##save this computation for use later
  if(ncol(H) != length(group)) {
    ## Contribution of other sites on the same choromosome
    pred_base <- pred_base + H[eval_idx, -group] %*% beta[-1,,drop=F][-group,,drop=F]
  }

  ## Compute test stat on original data
  if(f_stat) {
    true_results <- sapply(1:ncol(Y), function(k) {
      group.support <- which(beta[-1,k,drop=F][group,,drop=F]!=0)
      if(length(group.support)>0) {
        X0 <- haps_to_gen(H[eval_idx, group, drop=F][,group.support,drop=F])
        f.stat <- deviance_test_stats(haps_to_gen(pred_base)[,k], X0, Y[eval_y,k], family=family)
        return(f.stat)
      } else {
        return(0)
      }
    })
  } else {
    true_pred <- pred_base + H[eval_idx, group, drop=F] %*% beta[-1,,drop=F][group,,drop=F]
    if (family == "gaussian") {
      true_results <- col_cors(true_pred, Y[eval_y,,drop=F], genotypes = genotypes)
    } else if (family == "binomial") {
    true_results <- logistic_ll(true_pred, Y[eval_y,,drop=F], genotypes = genotypes)
    }
  }
  cat("done.\n")

  ## Select window
  window_size <- 200
  w.start <- group[1]
  w.end <- group[length(group)]
  p <- ncol(H)
  start_window <- max(1, w.start - window_size)
  end_window <- min(p, w.end + window_size)
  window <- start_window:end_window
  w.group <- which(window %in% group)

  worker.FUN <- function(i.rep) {
    ## Impute group within window
    H_impute <- t(vapply(eval_idx, function(i2) {
      digitaltwins::impute_variant_fast(H[i2,window], H_parent[anc[i2,1],window], H_parent[anc[i2,2],window],
                                      w.group, d[window], fb_result=fb_full[i2,], window_size=window_size)
    }, double(length(group))))
    ## Check whether impute_variant_fast failed (why would it fail?)
    if(nrow(H_impute) != length(eval_idx)) {
      return(rep(NA,ncol(Y)))
    }
    
    ##compute the CRT statistic on the replicate
    if(f_stat) {
      result <- sapply(1:ncol(Y), function(k) {
        group.support <- which(beta[-1,k,drop=F][group,,drop=F]!=0)
        if(length(group.support)>0) {
          X1 <- haps_to_gen(H_impute[,group.support,drop=F])
          f.stat.new <- deviance_test_stats(haps_to_gen(pred_base)[,k], X1, Y[eval_y,k], family=family)
          return(f.stat.new)
        } else {
          return(0)
        }
      })
    } else {
      pred <- pred_base + H_impute %*% beta[-1,,drop=F][group,,drop=F]
      if(family == "gaussian") {
        result <- col_cors(pred, Y[eval_y,,drop=F], genotypes)
      } else if (family == "binomial") {
        result <- logistic_ll(pred, Y[eval_y,,drop=F], genotypes)
      }
    }
    return(result)
  }

  if(parallel) {
    suppressMessages(library(doRNG))
    cat(sprintf("Computing CRT p-values over %d repetitions (multicore):\n", n_reps))
    results <- foreach(i = 1:n_reps, .combine = 'rbind') %dorng% { worker.FUN(i) }
  } else {
    cat(sprintf("Computing CRT p-values over %d repetitions (using single core):\n", n_reps))
    cat(sprintf("|%s|\n", paste(rep("-", 100), collapse = "")))
    cat("|")
    pbapply::pboptions(type = "txt", style = 1, char = "=", txt.width=100)
    results <- do.call("rbind", pbapply::pblapply(1:n_reps, worker.FUN))
  }

  ## Compute the p-values
  true_results <- matrix(true_results, nrow=n_reps, ncol=length(true_results), byrow=T)
  pvals <- colMeans(rbind(rep(TRUE,ncol(results)), true_results<=results), na.rm = TRUE)
  return(pvals)
}

#' DTT with independence, using a fixed linear statistic
#'
#' Returns p-values from the modified DTT using a linear feature importance measure. The p-values
#' are guaranteed to be independent.
#
# Arguments:
#' @param H an (2n x p) matrix of the subject haplotypes.
#'     It is assumed rows 1,2 belong to subject 1, rows 3,4 belong to subject 2, etc.
#' @param H_parent and (n_2 x p) matrix of the parental haplotypes
#' @param anc a (2n x 2) table of ancestries.
#'     anc[i, 1] and and[i, 2] give the rows of H_parent corresponding to the
#'     parents of row i of the haplotypes H
#' @param Y observed responses, a vector of length n or a matrix of dimension (n x k)
#' @param beta feature importance directions, fit on from an ind data set.
#'     A vector of length p+1 or a matrix of dimension ((p+1) x k).
#'     e.g. The output of get_beta_glment. The first coefficient is assumed to
#'     be an intercept.
#' @param d a vector of length p of genentic distances
#' @param adjust a vector of length n or a matrix of dimension (n x k),
#'     giving the contribution of the other chromosomes to the likelihood.
#'     I.e. adjust = X %*% beta_hat where X is all other chromosomes and beta is
#'     the fitted coefficients.
#' @param groups A list of of groups indices of the group to test. Each element should be 
#'     a continuous region, e.g. list(10:20, 21:30).
#' @param test_idx a set of indices (between 1 and 2n) used to compute the test statistics
#' @param family type of regression. Either "gaussian" or "binomial".
#'     If "guassian", correlation is used as a feature importance statistic. If "binomial",
#'    logistic log-likelihood is used as a feature importance statistic.
#' @param n_reps number of repitions of the CRT to carry out
#' @param genotypes Defaults to TRUE. If FALSE, all rows are assumed to be independent
#'     (i.e. no two haplotypes are from the same individual).
#' @param verbose if TRUE, prints various diagnostic messages to console
#' @param parallel requires doMC to be registered (default: FALSE)
#' @param recomb_thresh threshold of probability of recombination events to ignore. Lower
#'     values will have potentially higher power at increased computational cost.
#' @param f_stat whether or not to use the f-statistic for continuous responses.
#'     The f-statistic may be much slower for groups with many nonzero coefficients.
#
# Returns:
#' @return A matrix with k rows and length(groups) columns. 
#' Entry (i, j) is a p-value for group j using column i of Y as a response. 
#' @export
linear_crt_indep <- function(H, H_parent, anc, Y, beta, groups, d, adjust = NULL,
                       test_idx = 1:nrow(H), family = "gaussian", n_reps = 100,
                       genotypes = TRUE, verbose = FALSE, parallel = FALSE,
                       recomb_thresh = 0.001, f_stat = TRUE) {

  stopifnot( is.matrix(H) )
  stopifnot( is.matrix(H_parent) )
  stopifnot( is.matrix(anc) )

  ## Check dimensions of input
  stopifnot( ncol(H) == ncol(H_parent) )
  #stopifnot( (2*nrow(H)) == nrow(H_parent) )
  stopifnot( nrow(H) == nrow(anc) )
  stopifnot( ncol(anc) == 2 )
  #stopifnot( length(group) <= ncol(H) )
  stopifnot( length(d) == ncol(H) )
  stopifnot( length(test_idx) <= nrow(H) )

  if(is.null(dim(beta))) {
    beta <- as.matrix(beta, ncol = 1)
  }
  if(is.null(dim(Y))) {
    Y <- as.matrix(Y, nrow=length(Y), ncol = 1)
  }
  stopifnot( ncol(Y) == ncol(beta) )
  stopifnot( (2*nrow(Y)) == nrow(H) )
  stopifnot( nrow(beta) == (ncol(H)+1) )

  if(is.null(adjust)) {
    adjust <- matrix(0,nrow(Y),ncol(beta))
  }
  if(is.null(dim(adjust))) {
    adjust <- as.matrix(adjust,length(adjust),1)
  }
  stopifnot( nrow(adjust) == nrow(Y) )
  stopifnot( ncol(adjust) == ncol(Y) )

  #global null case
  if(length(groups) == 1 & length(groups[1] == ncol(H))) {
    return(1)
  }

  # sample the ancestry matrix, 1-indexed
  cat("Sampling the ancestry...")
  U = matrix(0, nrow = nrow(H), ncol = ncol(H))
  for(i in test_idx) {
    U[i, ] = digitaltwins:::sample_ancestry_(H[i, ], H_parent[anc[i, 1], ], H_parent[anc[i, 2], ], 
      p = ncol(H), epsilon = .001, lambda = .012, d) + 1
  }
  cat(" done.\n")
    
  full_pvals = c()
  cat(unique(groups))
  for(gid in unique(groups)) {
    if(gid == 0) {next} #flag for variables to ignore
    group = which(groups == gid)
    cat(sprintf("Starting group number: %d.\n", gid))
    cat(group)
    cat("\n")

    ## Ignore observations with no variability
    impute_idx <- which(U[, group[1]] != U[, group[length(group)]])
    cat(impute_idx)
    if(length(impute_idx) == 0) {
      cat("Warning, no recombination events detected.\n")
      full_pvals <- rbind(full_pvals, rep(1, ncol(beta)))
      next
    } else {
      cat(sprintf("Number of possible recombination events detected: %d.\n", length(impute_idx)))
    }

    ## Select ids to use in test statistic
    if(genotypes) {
      eval_idx <- get_genotype_ids(impute_idx)
      eval_y <- haps_id_to_gen(eval_idx)
    } else {
      eval_idx <- impute_idx
      eval_y <- eval_idx
    }

    cat("Computing stats... \n")
    ## Pre-compute some test stat information
    if(genotypes) {
      intercept <- matrix(beta[1,,drop=F], nrow=length(eval_idx), ncol=ncol(beta), byrow=TRUE) / 2
      adjust_mat <- adjust[rep(1:nrow(adjust),each=2),,drop=F][eval_idx,,drop=F] / 2
    } else {
      intercept <- matrix(beta[1,,drop=F], nrow=length(eval_idx), ncol=ncol(beta), byrow=TRUE)
      adjust_mat <- adjust[eval_idx,,drop=F]
    }
    pred_base <- intercept + adjust_mat ##save this computation for use later
    if(ncol(H) != length(group)) {
      ## Contribution of other sites on the choromosome
      pred_base <- pred_base + H[eval_idx, -group] %*% beta[-1,,drop=F][-group,,drop=F]
    }
    
    ## Compute test stat on original data
    if(f_stat) {
      true_results <- sapply(1:ncol(Y), function(k) {
        group.support <- which(beta[-1,k,drop=F][group,,drop=F]!=0)
        if(length(group.support)>0) {
          X0 <- haps_to_gen(H[eval_idx, group, drop=F][,group.support,drop=F])
          f.stat <- deviance_test_stats(haps_to_gen(pred_base)[,k], X0, Y[eval_y,k], family=family)
          return(f.stat)
        } else {
          return(0)
        }
      })
    } else {
      true_pred <- pred_base + H[eval_idx, group, drop=F] %*% beta[-1,,drop=F][group,,drop=F]
      if (family == "gaussian") {
        true_results <- col_cors(true_pred, Y[eval_y,,drop=F], genotypes = genotypes)
      } else if (family == "binomial") {
        true_results <- logistic_ll(true_pred, Y[eval_y,,drop=F], genotypes = genotypes)
      }
    }
    cat("done.\n")
         
    ## Select window
    start <- group[1]
    end <- group[length(group)]
    width = end - start + 1
    p <- ncol(H)

    worker.FUN <- function(i.rep) {
      ## Impute group within window
      H_impute  = matrix(0, ncol = length(group), nrow = length(eval_idx))
      for(idx in 1:length(eval_idx)) {
        i = eval_idx[idx]
        if(U[i, start] == U[i, end]) {
          H_impute[idx, ] = H[i, group]
        } else {

          # sample a recombination point weighted by genetic distance
          recomb_point = sample(1:(width - 1), 1, prob = d[I(start+1):end] / sum(d[I(start+1):end]))
          new_off = c(H_parent[anc[i, U[i, start]], group[1:recomb_point]],
                            H_parent[anc[i, U[i, end]], group[(recomb_point+1):width]])
          H_impute[idx, ] = new_off
        }
      }

      ##compute the CRT statistic on the replicate
      if(f_stat) {
        result <- sapply(1:ncol(Y), function(k) {
          group.support <- which(beta[-1,k,drop=F][group,,drop=F]!=0)
          if(length(group.support)>0) {
            X1 <- haps_to_gen(H_impute[,group.support,drop=F])
            f.stat.new <- deviance_test_stats(haps_to_gen(pred_base)[,k], X1, Y[eval_y,k], family=family)
            return(f.stat.new)
          } else {
            return(0)
          }
        })
      } else {
        pred <- pred_base + H_impute %*% beta[-1,,drop=F][group,,drop=F]
        if(family == "gaussian") {
          result <- col_cors(pred, Y[eval_y,,drop=F], genotypes)
        } else if (family == "binomial") {
          result <- logistic_ll(pred, Y[eval_y,,drop=F], genotypes)
        }
      }
      return(result)
    }

    if(parallel) {
      library(doRNG)
      cat(sprintf("Computing CRT p-values over %d repetitions (multicore):\n", n_reps))
      results <- foreach(i = 1:n_reps, .combine = 'rbind') %dorng% { worker.FUN(i) }
    } else {
      cat(sprintf("Computing CRT p-values over %d repetitions (using single core):\n", n_reps))
      results <- do.call("rbind", lapply(1:n_reps, worker.FUN))
    }

    ## Compute the p-values
    true_results <- matrix(true_results, nrow=n_reps, ncol=length(true_results), byrow=T)
    pvals <- colMeans(rbind(rep(TRUE,ncol(results)), true_results<=results), na.rm = TRUE)
    full_pvals <- rbind(full_pvals, pvals)
  }

  return(full_pvals)
}
