# Copyright (C) 2014 - 2018  Jack O. Wasey
#
# This file is part of icd.
#
# icd is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# icd is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with icd. If not, see <http:#www.gnu.org/licenses/>.

context("ICD-10 comorbidity calculations")
# since most of the code is common to ICD-9 and ICD-10, this doesn't need to be
# very extensive.

test_that("ICD-10 comorbidities from uranium", {
  expect_warning(comorbid(uranium_pathology, icd10_map_quan_elix), regexp = NA)
  expect_warning(comorbid(uranium_pathology, icd10_map_quan_deyo), regexp = NA)
  expect_warning(comorbid(uranium_pathology, icd10_map_elix), regexp = NA)
  expect_warning(comorbid(uranium_pathology, icd10_map_ahrq), regexp = NA)
})

context("icd10 comorbidity lookups")

test_that("ahrq comorbidities found for test data", {
  # test_two are all invalid codes
  test_two <- icd10_all_ahrq_one_pt
  test_two$icd10_code <- as.icd10cm(paste0(icd10_all_ahrq_one_pt$icd10_code, "J2"))
  test_three <- icd10_all_ahrq_one_pt
  test_three$icd10_code <- paste0(icd10_all_ahrq_one_pt$icd10_code,
                                  random_string(nrow(test_three), max_chars = 5))
  td <- named_list(icd10_all_ahrq_one_pt, test_two, test_three)
  for (test_name in names(td)) {
    expect_error(res <- comorbid(td[[test_name]], map = icd10_map_ahrq), regexp = NA, info = test_name)
    for (n in colnames(res))
      expect_true(res[, n], info = paste("method one comorbidity:", n, ", test: ", test_name))
    expect_error(res <- icd10_comorbid(td[[test_name]], map = icd10_map_ahrq), regexp = NA, info = test_name)
    for (n in colnames(res))
      expect_true(res[, n], info = paste("method two comorbidity:", n, ", test: ", test_name))
    expect_error(res <- icd10_comorbid_ahrq(td[[test_name]], hierarchy = FALSE, icd_name = "icd10_code"),
                 regexp = NA, info = test_name)
    for (n in colnames(res))
      expect_true(res[, n], info = paste("method three comorbidity:", n, ", test: ", test_name))

  }
})

test_that("ahrq cmb for many pts each with a diff cmb", {
  test_three <- test_two <- icd10_all_ahrq
  test_two$icd10_code <- paste0(icd10_all_ahrq$icd10_code, "J2")
  test_three$icd10_code <- paste0(icd10_all_ahrq$icd10_code,
                                  random_string(nrow(test_three), max_chars = 5))
  td <- list(
    icd10_all_ahrq = icd10_all_ahrq,
    test_two = test_two,
    test_three = test_three
  )

  # just checking for a trace matrix of TRUE, but this way lets me find error much more quickly:
  for (test_name in names(td)) {

    res <- comorbid(td[[test_name]], icd10_map_ahrq)
    for (n in colnames(res))
      expect_equal(sum(res[, n]), 1, info = paste("col method one comorbidity:", n, ", test: ", test_name))
    for (n in rownames(res))
      expect_equal(sum(res[n, ]), 1, info = paste("row method one comorbidity:", n, ", test: ", test_name))

    res <- icd10_comorbid(td[[test_name]], icd10_map_ahrq)
    for (n in colnames(res))
      expect_equal(sum(res[, n]), 1, info = paste("method two comorbidity:", n, ", test: ", test_name))
    for (n in rownames(res))
      expect_equal(sum(res[n, ]), 1, info = paste("row method two comorbidity:", n, ", test: ", test_name))
  }
})

test_that("cmb from a ICD-10, no infinite recursion with Elix", {
  # this code I5020 is in all the maps for ICD-10
  for (code in c("I5020")) {
    for (class_fun in c("as.character", "as.icd10", "as.icd10cm")) {
      for (map_fun in c("icd10_comorbid_elix", "icd10_comorbid_quan_elix",
                        "icd10_comorbid_quan_deyo", "icd10_comorbid_ahrq",
                        "comorbid_elix", "comorbid_quan_elix",
                        "comorbid_quan_deyo", "comorbid_ahrq")) {
        code_w_class <- do.call(class_fun, list(code))
        df <- data.frame(visit = 1, code = code_w_class, stringsAsFactors = FALSE)
        expect_error(res <- do.call(map_fun, list(x = df)), regex = NA, info = paste(class_fun, map_fun))
        expect_true(res[, "CHF"])
        expect_false("HTNcx" %in% names(res))
        df2 <- data.frame(visit = c(1, 2), code = c(code_w_class, code_w_class), stringsAsFactors = FALSE)
        expect_error(res <- do.call(map_fun, list(df2)), regex = NA, info = paste(class_fun, map_fun))
        expect_true(all(res[, "CHF"]))
        expect_false("HTNcx" %in% names(res))
      }
    }
  }
})

test_that("ICD-10 map reduction is sane", {
  uranium_short_codes <- as_char_or_levels(decimal_to_short.icd10(uranium_pathology$icd10))
  red_map <- simplify_map_lex(uranium_short_codes, icd10_map_ahrq)
  # the map should not have its original contents, but only those codes which
  # were in the pt data. Converse is not true, as some patient codes do not fall
  # into comorbidity.
  expect_true(all(unname(unlist(red_map)) %in% uranium_short_codes))
})

test_that("using reduction method for ICD-10", {
  if (!exists("icd10_comorbid_reduce"))
    skip("icd10_comorbid_reduce not available")
  res <- icd10_comorbid(uranium_pathology, map = icd10_map_ahrq, icd10_comorbid_fun = icd10_comorbid_reduce)
  expect_equal(ncol(res), 30)

  expect_identical(
    icd10_comorbid_parent_search_use_cpp(uranium_pathology, icd10_map_ahrq,
                                         visit_name = "case", icd_name = "icd10",
                                         short_code = FALSE, short_map = TRUE, return_df = FALSE),
    icd10_comorbid_reduce(uranium_pathology, icd10_map_ahrq,
                          visit_name = "case", icd_name = "icd10",
                          short_code = FALSE, short_map = TRUE, return_df = FALSE))

})
