# conditional_tdt
#' The DTT with a TDT-like test statistic
#'
#' Runs the local Digital Twin Test using the TDT test statistc. This test
#' determines if there is a causal SNP somewhere in a specified windown, after
#' accounting for all genetic variant outside the window.
#' 
#' This function is not intented to be called repeatedly for different variants. 
#' Significant computational speedups would be possible in that case.
#'
# Arguments:
#' @param sample A dataframe in .sample format. \code{sample$father} gives the index of
#'      the father in the \code{haps_matrix}, and the \code{sample$mother} gives the index of the 
#'      the mother in the \code{haps_matrix}.
#' @param haps_matrix A n x p matrix of 0-1. Each row represents one haplotypes. 
#'      To convert from indices to haplotype indices, an entry i in \code{sample$father}
#'      corresponds to the two rows 2*i - 1 and 2*i in \code{haps_matrix}.
#' @param snp_info A dataframe with 5 columns: chromosome, snp_name, 
#'      position, variant_1, and variant_2 (column naming is ignored).
#' @param group A vector of entries specifying the group to test. Must be a vector 
#'      of the form j:k for some j > 0 and k < p + 1.
#' @param site The index of a column of \code{haps_matrix} to test. 
#'      This is analagous to the variant tested by the TDT.
#' @param gen_map Genetic map: a dataframe contining at least columns "pposition" and "gposition".
#' @param n_reps Number of repetions of the randomization test.
#' @param exact Whether or not to include the finite-sample correction.
#
# Returns:
#' @return A p-value.
#' @export
conditional_tdt = function(sample, haps_matrix, snp_info, group, site, gen_map, 
		n_reps = 500, exact = FALSE) {
    cat("Resolve maternal vs paternal haplotype identities... ")
	#resolve maternal vs paternal haplotype identity
	alignments = get_phased_ancestry(sample, haps_matrix)
	#select full trios
	alignments = dplyr::filter(alignments, ancestor_2 != "")
    cat("done.\n")

	#find test and parent indices
    test_idx_gen = match(as.character(alignments$subject), as.character(sample$ID_2))
    test_idx = rep(0, 2*length(test_idx_gen))
    test_idx[(1:length(test_idx_gen))*2 - 1] = 2*test_idx_gen - 1
    test_idx[(1:length(test_idx_gen))*2 ] = 2*test_idx_gen
    anc_1 = match(as.character(alignments$ancestor_1), as.character(sample$ID_2))
    anc_2 = match(as.character(alignments$ancestor_2), as.character(sample$ID_2))
    anc = matrix(0, nrow = length(test_idx), ncol = 2)
    anc[1:(nrow(anc)/2)*2 - 1, ] = cbind(2 * anc_1 - 1, 2 * anc_1)
    anc[1:(nrow(anc)/2)*2, ] = cbind(2 * anc_2 - 1, 2 * anc_2)
    stopifnot(length(test_idx) == nrow(anc))

	#extract information from the genetic map
	colnames(snp_info) = c("chromosome", "snp_name", "position", "V1", "V2")
	if(!is.null(gen_map)) {
		temp = dplyr::left_join(snp_info, gen_map, by = c("position" = "pposition"))
		#linearly fill missing genetic position of missing entries
		g_pos = temp$gposition
        n_snp = length(g_pos)
		g_pos[1] = 0
		g_pos[n_snp] = max(c(0, g_pos), na.rm = TRUE)
        observed = which(!is.na(g_pos))
		for(j in 2:(n_snp - 1)) {
			if(is.na(g_pos[j])) {
				next_observed = min(observed[observed > j])
				gap = next_observed - j + 1
				g_pos[j] = 1 / gap * (g_pos[next_observed] - g_pos[j - 1]) + g_pos[j - 1]
			}
		}
	}
	d = g_pos[2:n_snp] - g_pos[1:(n_snp - 1)]
	cat("Fraction of sites in genetic map: ", length(observed) / length(g_pos), "\n")

	#forward backward pass to determine ancestry of the observed haplotype
    cat("Pre-computing forward-backward results... ")
    fb_results <- get_fb_results(haps_matrix[test_idx, ], haps_matrix, anc = anc, 
    							 group = group, d = d,
                                 lambda = .01, epsilon = .001, window_size = 200)
    cat("done.\n")

    #extract probabilities for the conditional TDT
    pos_in_group = (g_pos[site] - g_pos[min(group)]) / (g_pos[max(group)] - g_pos[min(group)])
    #conditional probability of inheriting from first haplotype of ancestor
    p1 = fb_results[, 1] * (1 - fb_results[, 2]) +
    	fb_results[, 1] * fb_results[, 2] * (1 - pos_in_group) + 
    	(1 - fb_results[, 1]) * fb_results[, 3] * pos_in_group

    #carry out the conditional TDT
    cat("Carrying out the randomizations... ")
    t_obs = sum(haps_matrix[test_idx, site]) #value on observed data
    t_sim = c() #value on simulated replicates
    for(k in 1:n_reps) {
    	u = runif(length(test_idx))
    	t_temp = sum(haps_matrix[anc[u < p1, 1], site]) + 
    		sum(haps_matrix[anc[u > p1, 2], site])
    	t_sim = c(t_sim, t_temp)
    }
    cat("done.\n")
    pv = 2 * min(mean(t_obs > t_sim), mean(t_obs < t_sim))
    if(exact) {pv = pv + 2 / n_reps}

    return(pv)
}