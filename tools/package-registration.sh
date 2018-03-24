#!/bin/bash
Rscript -e 'tools::package_native_routine_registration_skeleton(".", "src/registration.c", character_only = FALSE)'
# try to drop unwanted registrations, but this led to other problems, as Rcpp::compileAttributes still generated code which needed registration
# sed -i '' '/_alt_/d' src/registration.c

