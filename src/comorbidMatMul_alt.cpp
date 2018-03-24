// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r, cpp)]]
// #ifdef ICD_EIGEN
#include <Rcpp.h>
#include <RcppEigen.h> // also add LinkingTo element in DESCRIPTION to enable
#include <Eigen/SparseCore>
#include "comorbidCommon.h"
#include "comorbidSetup.h"
#include <algorithm>                   // for binary_search, copy
#include <vector>                      // for vector, vector<>::const_iterator
#include "icd_types.h"                 // for ComorbidOut, VecVecInt, VecVec...
#include "local.h"                     // for ICD_OPENMP
#include "config.h"                     // for valgrind, CXX11 etc
#include "util.h"                     // for debug_parallel
extern "C" {
#include "cutil.h"                              // for getRListOrDfElement
}

using namespace Rcpp;

// use row-major sparse matrix - row major because easier to insert into Eigen
// sparse matrix, and we discover comorbidities one patient at a time, i.e. row
// major

// using the typedef confuses Rcpp
//typedef Eigen::SparseMatrix<char, Eigen::RowMajor> SparseOut; // bool, char or int?
// https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
typedef Eigen::Triplet<char> Triplet;
typedef Eigen::SparseMatrix<char, Eigen::RowMajor> PtsSparse;

//' @title prototype to do entire comorbidity calculation as a matrix multiplication
//' @description
//' The problem is that the matrices could be huge: the patient-icd matrix would
//' be millions of patient rows, and ~15000 columns for all AHRQ comorbidities.
//' @details
//' Several ways of reducing the problem: firstly, as with existing code, we can
//' drop any ICD codes from the map which are not in the patient data. With many
//' patients, this will be less effective as the long tail becomes apparent.
//' However, with the (small) Vermont data, we see ~15,000 codes being reduced to
//' 339.
//' @section Sparse matrices
//' Using sparse matrices is another solution. Building
//' the initial matrix may become a significant part of the calculation, but once
//' done, the solution could be a simple matrix multiplication, which is
//' potentially highly optimized (Eigen, BLAS, GPU, etc.)
//' @section Eigen
//' Eigen has parallel (non-GPU) optimized sparse row-major *
//' dense matrix. Patients-ICD matrix must be the row-major sparse one, so the
//' dense matrix is then the comorbidity map
//' @examples
//' # show how many discrete ICD codes there are in the AHRQ map, before reducing
//' # to the number which actually appear in a group of patient visits
//' sapply(icd::icd9_map_ahrq, length) %>% sum
//' icd_comorbid_ahrq(vermont_dx %>% icd_wide_to_long, comorbid_fun = icd:::icd9ComorbidShortCpp)
//' \dontrun{
//' # remove _alt line in .Rbuildignore, then these will be available. Also, re-enable [[Rcpp::depends(RcppEigen)]]
//' icd_comorbid_ahrq(vermont_dx %>% icd_wide_to_long, comorbid_fun = icd:::icd9ComorbidMatMul)
//' icd_comorbid_ahrq(vermont_dx %>% icd_wide_to_long, comorbid_fun = icd:::icd9ComorbidSparseOmp)
//' }
//' @keywords internal
// [[Rcpp::export]]
SEXP icd9Comorbid_alt_MatMul(const SEXP& icd9df, const Rcpp::List& icd9Mapping,
                             const std::string visitId, const std::string icd9Field,
                             const int threads = 8, const int chunk_size = 256,
                             const int omp_chunk_size = 1, bool aggregate = true) {

  // find eventual size of map matrix:
  size_t map_rows = 0;
  for (List::const_iterator li = icd9Mapping.begin(); li != icd9Mapping.end(); ++li) {
    IntegerVector v = *li;
    map_rows += v.size();
  }

  size_t row = 0;
  // make an integer matrix for the map. not sparse. No boolean option, I don't think.
  Eigen::MatrixXi map(map_rows, icd9Mapping.length());
  for (List::const_iterator li = icd9Mapping.begin(); li != icd9Mapping.end(); ++li) {
    IntegerVector v = *li;
    for (IntegerVector::iterator vi = v.begin(); vi != v.end(); ++vi, ++row) {
      map.coeffRef(row, std::distance(v.begin(), vi)) = true;
    }
  }
  Rcpp::Rcout << "first cell of map should be 1: it is " << map(0, 0) << "\n";

  // now build the patient:icd matrix... can probably re-use and simplify the
  // function 'buildVisitCodesVec'
  PtsSparse pts; // rows = unique patients, cols = unique ICD codes

  // Eigen::MatrixXi map(map_rows,
  //icd9Mapping.length());

  return List();
}
// #endif // ICD_EIGEN
