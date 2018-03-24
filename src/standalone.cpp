#ifdef INSIDE

// don't include my source files: this is like an external application so link against my Rcpp generated headers.
#include "icd.h"
#include <RcppCommon.h>
#include <Rcpp.h>
#include <RInside.h>
#include <string>
#include <iostream>


//int main(int argc, char *argv[]) {
int main() {
  using namespace Rcpp;

	RInside R(argc, argv);          // create an embedded R instance

	const Rcpp::CharacterVector testmaj = Rcpp::CharacterVector::create("E10", "V1", "5");
	//CV res = icd9GetMajor(testmaj, false);
	Rcpp::CharacterVector res;
	res = icd::icd9_add_leading_zeroes_major(testmaj);

	Rcpp::Rcout << "res = " << res << "\n";

	exit(0);
}

#endif // #INSIDE
