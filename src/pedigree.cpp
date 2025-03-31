// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Extend (or construct) the pedigree matrix for a generation
//'
//' Constructs a pedigree matrix for the provided generation of size (L * n)
//' where L is the pedigree depth and n is the population size. If the current
//' pedigree depth is L, the oldest layer is discarded. If the current pedigree
//' is an empty matrix, initialises the pedigree.
//'
//' @param pedigree_current The pedigree of the current generation to extend
//' @param matching Integer vector of length n, where matching[i] specifies the
//'   index of the mate of individual i in the population
//' @param pedigree_depth
//'
//' @return An extended pedigree for the current generation
//'
//' @noRd
arma::mat construct_pedigree(
  arma::mat pedigree_current,
  const arma::uvec matching,
  const unsigned int pedigree_depth = 3
) {
  return arma::mat(pedigree_current.n_rows, 2);
}

//' Compute the pedigree distance between individuals
//'
//' Computes the distance between two individuals in the population pedigree
//' graph given their index identifiers and the pedigree using the LCA algorithm
//' for DAGs. Returns the maximum integer value if a common ancestor is not
//' found.
//'
//' @param pedigree The population pedigree
//' @param ind1_id Index of the first individual
//' @param ind2_id Index of the second individual
//'
//' @return The pedigree distance between ind1 and ind2, or the maximum integer
//'   value if no common ancestor is found.
//'
//' @noRd
int compute_pedigree_distance(
  arma::mat pedigree,
  int ind1_id,
  int ind2_id
) {
  return std::numeric_limits<int>::max();
}
