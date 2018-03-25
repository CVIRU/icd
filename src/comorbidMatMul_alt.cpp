// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r, cpp)]]
#include "config.h"                     // for valgrind, CXX11 etc
#include "local.h"                     // for ICD_OPENMP
#ifdef ICD_EIGEN // rest of file
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
typedef int SparseValue;
typedef Eigen::Triplet<SparseValue> Triplet;
typedef Eigen::SparseMatrix<SparseValue, Eigen::RowMajor> PtsSparse;
typedef Eigen::MatrixXi DenseMap;

// alternate version which builds a sparse matrix, row-major, which is good for
// LHS of multrix multiplication in Eigen
void buildVisitCodesVecSparse(const SEXP& icd9df,
                              const std::string& visitId,
                              const std::string& icd9Field,
                              PtsSparse& sparse_db,
                              VecStr& visitIds, // will have to get this from sparse matrix at end, but needed?
                              const bool aggregate = true // remove or ignore
) {

  SEXP icds = PROTECT(getRListOrDfElement(icd9df, icd9Field.c_str())); // this is a factor
  SEXP vsexp = PROTECT(getRListOrDfElement(icd9df, visitId.c_str()));
  IntegerVector icd9dfFactor = as<IntegerVector>(icds);
  CV factorLevels = icd9dfFactor.attr("levels");
  R_xlen_t numUniqueCodes = factorLevels.length();

  const char* lastVisitId = "JJ94967295JJ"; // random

  int vlen = Rf_length(icds); // or vsexp
  // make an unordered set for quick check for duplicates while building list of unique visit ids
  std::unordered_set<std::string> visit_lookup;
  visit_lookup.reserve(vlen);
  // also maintain list of (ordered as first encountered) visit ids
  visitIds.resize(vlen); // resize and trim at end, as alternative to reserve
  int n;
  VisLk vis_lookup;
  VecVecIntSz vcdb_max_idx = -1; // we increment immediately to zero as first index
  VecVecIntSz vcdb_new_idx;
  VecVecIntSz vcdb_last_idx;

  sparse_db.resize(vlen, numUniqueCodes);
  sparse_db.reserve(vlen); // but memory commitment is known and limited.

  std::vector<Triplet> visTriplets;
  visTriplets.reserve(vlen * 30); // overestimate codes per patient to avoid resizing while filling

  for (int i = 0; i != vlen; ++i) {
    const char* visit = CHAR(STRING_ELT(vsexp, i));
    n = INTEGER(icds)[i]; // ICD codes are in a factor, so get the integer index

    if (lastVisitId != visit) {
      // assume new visitId unless aggregating
      vcdb_new_idx = vcdb_max_idx + 1;
      if (aggregate) { // only use map if aggregating
        VisLk::iterator found = vis_lookup.find(visit); // did we see this visit already?
        if (found != vis_lookup.end()) {
          // we saved the index in the map, so use that to insert a triplet:
          visTriplets.push_back(Triplet(found->second, n, true));
          continue; // and continue with next row
        } else { // otherwise we found a new visitId, so add it to our lookup table
          vis_lookup.insert(
            std::make_pair(visit, vcdb_new_idx)); // new visit, with associated position in vcdb
        }
      } // end if aggregate
      // we didn't have an existing visitId, or we are just assuming visitIds are ordered (not aggregating)
      visTriplets.push_back(Triplet(vcdb_new_idx, n, true));
      visitIds[vcdb_new_idx] = visit; // keep list of visitIds in order encountered.
      lastVisitId = visit;
      vcdb_last_idx = vcdb_new_idx;
      ++vcdb_max_idx;
    } else { // last visitId was the same as the current one, so we skip all the logic
      visTriplets.push_back(Triplet(vcdb_last_idx, n, true));
    }

  } // end loop through all visit-code input data
  UNPROTECT(2);

  // sparse_db and visitIds are updated
  sparse_db.setFromTriplets(visTriplets.begin(), visTriplets.end());
  visitIds.resize(vcdb_max_idx + 1); // we over-sized (not just over-reserved) so now we trim.
}

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
//' \url{https://eigen.tuxfamily.org/dox/TopicMultiThreading.html}
//' @examples
//' # show how many discrete ICD codes there are in the AHRQ map, before reducing
//' # to the number which actually appear in a group of patient visits
//' sapply(icd::icd9_map_ahrq, length) %>% sum
//' icd_comorbid_ahrq(vermont_dx %>% icd_wide_to_long, comorbid_fun = icd:::icd9ComorbidShortCpp)
//' \dontrun{
//' # remove _alt line in .Rbuildignore, then these will be available. Also, re-enable [[Rcpp::depends(RcppEigen)]]
//' icd_comorbid_ahrq(vermont_dx %>% icd_wide_to_long, comorbid_fun = icd:::icd9Comorbid_alt_MatMul)
//' icd_comorbid_ahrq(vermont_dx %>% icd_wide_to_long, comorbid_fun = icd:::icd9Comorbid_alt_SparseOmp)
//' }
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix icd9Comorbid_alt_MatMul(const SEXP& icd9df, const Rcpp::List& icd9Mapping,
                             const std::string visitId, const std::string icd9Field,
                             const int threads = 8, const int chunk_size = 256,
                             const int omp_chunk_size = 1, bool aggregate = true) {
  valgrindCallgrindStart(true);
  VecStr out_row_names; // size is reserved in buildVisitCodesVec
  // find eventual size of map matrix:
  size_t map_rows = 0;
  for (List::const_iterator li = icd9Mapping.begin(); li != icd9Mapping.end(); ++li) {
    IntegerVector v = *li;
    map_rows += v.size();
  }

  size_t row = 0;
  // make an integer matrix for the map. not sparse. No boolean option, I don't think.
  DenseMap map(map_rows, icd9Mapping.length());
  for (List::const_iterator li = icd9Mapping.begin(); li != icd9Mapping.end(); ++li) {
    IntegerVector v = *li;
    for (IntegerVector::iterator vi = v.begin(); vi != v.end(); ++vi, ++row) {
      map.coeffRef(row, std::distance(v.begin(), vi)) = true;
    }
  }
  Rcpp::Rcout << "first cell of map should be 1: it is " << map(0, 0) << std::endl;

  // now build the patient:icd matrix... can probably re-use and simplify the
  PtsSparse sparse_db; // reservation and sizing done next
  buildVisitCodesVecSparse(icd9df, visitId, icd9Field, sparse_db, out_row_names, aggregate);
  Rcpp::Rcout << " built the sparse matrix " << std::endl;
  DenseMap result = sparse_db * map;
  Rcpp::Rcout << " done matrix multiplication " << std::endl;
  Rcpp::IntegerMatrix mat_out = Rcpp::wrap(result);
  Rcpp::Rcout << " cast to integer matrix " << std::endl;
  mat_out.attr("dimnames") = Rcpp::List::create(icd9Mapping.names(), out_row_names);
  Rcpp::Rcout << "dimension names set" << std::endl;

  valgrindCallgrindStop();
  return mat_out;
}

#endif // ICD_EIGEN
