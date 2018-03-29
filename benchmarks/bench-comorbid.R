# generate a bunch of random patients and get comorbidities, intended for valgrind usage
library(icd)
library(microbenchmark)

pts <- icd:::generate_random_pts(5e6)
message("benchmark starting!")
mb <- microbenchmark(
  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_MatMul),
  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_Taskloop),
  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_Taskloop2),
  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9ComorbidShortCpp),
  times = 5L)

#  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_SparseOmp),

print(mb)

