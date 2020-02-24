# get_phased ancestry
#' Indentify maternal vs paternal haplotypes
#'
#' Identifies which of a subjects two haplotyes
#' is the maternal haplotype and which is the paternal haplotype.
#
# Arguments:
#' @param sample A dataframe in .sample format. \code{sample$father} gives the index of
#'      the father in the \code{haps_matrix}, and the \code{sample$mother} gives the index of the 
#'      the mother in the \code{haps_matrix}.
#' @param haps_matrix A n x p matrix of 0-1. Each row represents one haplotypes. 
#'      To convert from indices to haplotype indices, an entry i in \code{sample$father}
#'      corresponds to the two rows 2*i - 1 and 2*i in \code{haps_matrix}.
#
# Returns:
#' @return A data frame with the following elements:
#' @return subject: the row ID of the subject. An entry of i corresponds to
#'  rows 2*i -1 and 2*i of \code{haps_matrix}.
#' @return chromosome: the chromosome number 1-22.
#' @return ancestor_1: the ID of the ancestor corresponding to the first haplotype
#'      of the offspring in \code{haps_matrix}.
#' @return ancestor_2: the ID of the ancestor corresponding to the second haplotype
#'      of the offspring in \code{haps_matrix}.
get_phased_ancestry = function(sample, haps_matrix) {
    alignments = data.frame()

    for(i in 1:nrow(sample)) {
        if(sample$father[i] != 0 & sample$mother[i] != 0) {
            f_row = as.integer(rownames(sample)[as.character(sample$ID_2) == 
                                                  as.character(sample$father[i])]) - 1
            m_row = as.integer(rownames(sample)[as.character(sample$ID_2) == 
                                                  as.character(sample$mother[i])]) - 1
            
            hap_row = 2*i - 1
            f_hap_row = 2*f_row - 1
            m_hap_row = 2*m_row - 1
            
            #offspring's 1st haplotype
            f_match1 = mean(haps_matrix[hap_row, ] == haps_matrix[f_hap_row, ] | 
                            haps_matrix[hap_row, ] == haps_matrix[f_hap_row + 1, ])
            m_match1 = mean(haps_matrix[hap_row, ] == haps_matrix[m_hap_row, ] | 
                            haps_matrix[hap_row, ] == haps_matrix[m_hap_row + 1, ])
            #offspring's 2nd haplotype
            f_match2 = mean(haps_matrix[hap_row + 1, ] == haps_matrix[f_hap_row, ] | 
                            haps_matrix[hap_row + 1, ] == haps_matrix[f_hap_row + 1, ])
            m_match2 = mean(haps_matrix[hap_row + 1, ] == haps_matrix[m_hap_row, ] | 
                            haps_matrix[hap_row + 1, ] == haps_matrix[m_hap_row + 1, ])
            
            #extract the order of the parental haplotypes
            if(length(f_row) == 0 & length(m_row) == 0) {
                first_hap = ""
                second_hap = ""
            } else if (length(f_row) == 0) {
                first_hap = as.character(sample$mother[i])
                second_hap = ""       
            } else if (length(m_row) == 0) {
                first_hap = as.character(sample$father[i])
                second_hap = ""       
            } else if(f_match1 > m_match1 & f_match2 < m_match2 & 
               max(f_match1, m_match1) > .99 & max(f_match2, m_match2) > .99) {
                first_hap = as.character(sample$father[i])
                second_hap = as.character(sample$mother[i])
            } else if(f_match1 < m_match1 & f_match2 > m_match2 & 
               max(f_match1, m_match1) > .99 & max(f_match2, m_match2) > .99) {
                first_hap = as.character(sample$mother[i])
                second_hap = as.character(sample$father[i])
            } else {
                print(paste0("High number of mendelian inconsistencies for subject: ", 
                             sample$ID_2[i], ". Ancestry not resolved."))
                first_hap = ""
                second_hap = ""
            }
                
            alignments = rbind(alignments, c(as.character(sample$ID_2[i]), 22, first_hap, second_hap), 
                              stringsAsFactors = FALSE)
            
        }
    }
        
    colnames(alignments) = c("subject", "chromosome", "ancestor_1", "ancestor_2")

    alignments
}