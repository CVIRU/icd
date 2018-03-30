# generate a bunch of random patients and get comorbidities, intended for valgrind usage
library(icd)
library(microbenchmark)
library(R.cache)
ten_million_random_pts <- R.cache::loadCache(key = list("ten_million_random_pts"), suffix = "icd.Rcache")
if (is.null(ten_million_random_pts)) {
    ten_million_random_pts <- icd:::generate_random_pts(1e7)
    R.cache::saveCache(ten_million_random_pts, key = list("ten_million_random_pts"))
}
# the following is a mix of repeated vermont patients and random patients,
# totally about 114,000,000 rows, which is about 3.5GB on disk in R.cache and
# 9GB in RAM.
huge_mixed_pts <- R.cache::loadCache(key = list("huge_mixed_pts"), suffix = "icd.Rcache")
if (is.null(huge_mixed_pts)) {
  vt <- icd_wide_to_long(vermont_dx)[c("visit_id", "icd_code")]
  vts <- mefa:::rep.data.frame(vt, 10000)
  rnd <- rnd[c("visit_id", "code")]
  names(rnd) <- names(vts)
  huge_mixed_pts <- rbind(rnd, vts)
  R.cache::saveCache(huge_mixed_pts, key = list("huge_mixed_pts"))
}

message("benchmark starting!")
mb <- microbenchmark(
  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_MatMul),
  #  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_Taskloop),
#  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_Taskloop2),
  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9ComorbidShortCpp),
  times = 25L)

#  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_SparseOmp),

print(mb)

