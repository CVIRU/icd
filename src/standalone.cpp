//#ifdef ICD_STANDALONE

// don't include my source files: this is like an external application so link against my Rcpp generated headers.
#include <RInside.h>
#include "icd.h"
#include <string>
#include <iostream>

int main(int argc, char *argv[]) {

	RInside R(argc, argv);          // create an embedded R instance

	const Rcpp::CharacterVector testmaj = Rcpp::CharacterVector::create("E10", "V1", "5");
	//CV res = icd9GetMajor(testmaj, false);
	Rcpp::CharacterVector res;
	res = icd::icd9_add_leading_zeroes_major(testmaj);

	Rcpp::Rcout << "res = " << res << "\n";

	return 0; // or exit(0) ?
}

//#endif // #STANDALONE
