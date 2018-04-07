# generate a bunch of random patients and get comorbidities, intended for valgrind usage
library(icd)
library(microbenchmark)
library(R.cache)
library(profvis)
ten_million_random_pts <- R.cache::loadCache(key = list("ten_million_random_pts"), suffix = "icd.Rcache")
if (is.null(ten_million_random_pts)) {
  ten_million_random_pts <- icd:::generate_random_pts(1e7)
  R.cache::saveCache(ten_million_random_pts, key = list("ten_million_random_pts"), suffix = "icd.Rcache")
}
# the following is a mix of repeated vermont patients and random patients,
# totally about 114,000,000 rows, which is about 3.5GB on disk in R.cache and
# 9GB in RAM.
huge_mixed_pts <- R.cache::loadCache(key = list("huge_mixed_pts"), suffix = "icd.Rcache")
if (is.null(huge_mixed_pts)) {
  vt <- icd_wide_to_long(vermont_dx)[c("visit_id", "icd_code")]
  vts <- mefa:::rep.data.frame(vt, 10000)
  rnd <- ten_million_random_pts[c("visit_id", "code")]
  names(rnd) <- names(vts)
  huge_mixed_pts <- rbind(rnd, vts)
  huge_mixed_pts$visit_id <- icd:::as_char_no_warn(huge_mixed_pts$visit_id)
  rm(list = c("rnd", "vts", "vt"))
  R.cache::saveCache(huge_mixed_pts, key = list("huge_mixed_pts"), suffix = "icd.Rcache")
}

message("profiling!")
profvis(icd_comorbid_ahrq(ten_million_random_pts, preclean = FALSE))

profvis::profvis(
  icd::icd_comorbid_ahrq(huge_mixed_pts, preclean = FALSE)
)

message("benchmark starting!")
mb <- microbenchmark(
  icd_comorbid(huge_mixed_pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_MatMul),
  #  icd_comorbid(huge_mixed_pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_Taskloop),
  #  icd_comorbid(huge_mixed_pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_Taskloop2),
  icd_comorbid(huge_mixed_pts, icd9_map_ahrq, comorbid_fun = icd:::icd9ComorbidShortCpp),
  times = 25L)

#  icd_comorbid(pts, icd9_map_ahrq, comorbid_fun = icd:::icd9Comorbid_alt_SparseOmp),

print(mb)

big_icd10 <- do.call("rbind", replicate(250, uranium_pathology, simplify = FALSE))

# look at all icd-10 options at once
microbenchmark::microbenchmark(
  icd:::icd10_comorbid_reduce(
    uranium_pathology, icd10_map_ahrq,
    visit_name = "case", icd_name = "icd10",
    short_code = FALSE, short_map = TRUE, return_df = FALSE),
  icd:::icd10_comorbid_parent_search_use_cpp(
    uranium_pathology, icd10_map_ahrq,
    visit_name = "case", icd_name = "icd10",
    short_code = FALSE, short_map = TRUE, return_df = FALSE),
  check = icd:::all_identical, times = 30)
# The following are all much slower by 1-2 orders of magnitude

# icd:::icd10_comorbid_parent_search_str(
#   uranium_pathology, icd10_map_ahrq,
#   visit_name = "case", icd_name = "icd10",
#   short_code = FALSE, short_map = TRUE, return_df = FALSE),
# icd:::icd10_comorbid_parent_search_orig(
#   uranium_pathology, icd10_map_ahrq,
#   visit_name = "case", icd_name = "icd10",
#   short_code = FALSE, short_map = TRUE, return_df = FALSE),
# icd:::icd10_comorbid_parent_search_all(
#   uranium_pathology, icd10_map_ahrq,
#   visit_name = "case", icd_name = "icd10",
#   short_code = FALSE, short_map = TRUE, return_df = FALSE),
# icd:::icd10_comorbid_parent_search_no_loop(
#   uranium_pathology, icd10_map_ahrq,
#   visit_name = "case", icd_name = "icd10",
#   short_code = FALSE, short_map = TRUE, return_df = FALSE),

# Found tight confidence intervals on reduce method being 1-2x as fast as cpp parent
# search for small (uranium) data

microbenchmark::microbenchmark(
  icd:::icd10_comorbid_reduce(
    big_icd10, icd10_map_ahrq,
    visit_name = "case", icd_name = "icd10",
    short_code = FALSE, short_map = TRUE, return_df = FALSE),
  icd:::icd10_comorbid_parent_search_use_cpp(
    big_icd10, icd10_map_ahrq,
    visit_name = "case", icd_name = "icd10",
    short_code = FALSE, short_map = TRUE, return_df = FALSE),
  identical = icd:::my_check, times = 10)

# reduce method 1.5-2 orders of magnitude faster!
