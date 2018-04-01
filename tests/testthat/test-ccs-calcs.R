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

context("ahrq ccs calculations")

test_that("icd9 CCS map is valid", {
  expect_true(icd:::icd_is_valid.icd_comorbidity_map(icd9_map_single_ccs, short_code = TRUE))
})

test_that("one code from each single level", {
  first_from_each <-
    vapply(icd9_map_single_ccs, function(y) as.character(y[[1]]), character(1), USE.NAMES = FALSE)
  # drop the (maybe unnecessary empty first group)
  first_from_each <- first_from_each[first_from_each != ""]
  test_all_ccs_df <- data.frame(
    visit_id = rep("z", length(first_from_each)),
    icd9 = first_from_each
  )
  res <- icd9_comorbid_ccs(test_all_ccs_df,  visit_name = "visit_id", icd_name = "icd9")
  res2 <- icd9_comorbid_ccs(test_all_ccs_df,  visit_name = "visit_id", icd_name = "icd9", comorbid_fun = icd:::icd9ComorbidShortCpp)
  expect_identical(res, res2)
  # should be one for each (do this way to ignore the empty first group)
  expect_equal(sum(res), length(first_from_each))

  # same but with strings instead of factors
  test_all_ccs_df <- data.frame(
    visit_id = rep("z", length(first_from_each)),
    icd9 = first_from_each,
    stringsAsFactors = FALSE
  )
  res <- icd9_comorbid_ccs(test_all_ccs_df,  visit_name = "visit_id", icd_name = "icd9")
  res2 <- icd9_comorbid_ccs(test_all_ccs_df,  visit_name = "visit_id", icd_name = "icd9", comorbid_fun = icd:::icd9ComorbidShortCpp)
  expect_identical(res, res2)
  # should be one for each (do this way to ignore the empty first group)
  expect_equal(sum(res), length(first_from_each))
})

test_that("one code from each single level backwards", {
  first_from_each <- rev(
    vapply(icd9_map_single_ccs, function(y) as.character(y[[1]]), character(1), USE.NAMES = FALSE))
  # drop the (maybe unnecessary empty first group)
  first_from_each <- first_from_each[first_from_each != ""]
  test_all_ccs_df <- data.frame(
    visit_id = rep("z", length(first_from_each)),
    icd9 = first_from_each
  )
  res <- icd9_comorbid_ccs(test_all_ccs_df,  visit_name = "visit_id", icd_name = "icd9")
  res2 <- icd9_comorbid_ccs(test_all_ccs_df,  visit_name = "visit_id", icd_name = "icd9", comorbid_fun = icd:::icd9ComorbidShortCpp)
  expect_identical(res, res2)
  # should be one for each (do this way to ignore the empty first group)
  expect_equal(sum(res), length(first_from_each))
})

test_that("one code from each single level backwards with disordered visits", {
  first_from_each <- rev(
    vapply(icd9_map_single_ccs, function(y) as.character(y[[1]]), character(1), USE.NAMES = FALSE))
  # drop the (maybe unnecessary empty first group)
  first_from_each <- first_from_each[first_from_each != ""]
  set.seed(1441)
  rnd_ccs_df <- data.frame(
    visit_id = sample(c("j", "b", "k"), size = length(first_from_each), replace = TRUE),
    icd9 = first_from_each
  )
  res <- icd9_comorbid_ccs(rnd_ccs_df,  visit_name = "visit_id", icd_name = "icd9")
  res2 <- icd9_comorbid_ccs(rnd_ccs_df,  visit_name = "visit_id", icd_name = "icd9", comorbid_fun = icd:::icd9ComorbidShortCpp)
  expect_identical(res, res2)
  # should be one for each (do this way to ignore the empty first group)
  expect_equal(sum(res), length(first_from_each))
})

test_that("smaller test case based on random input", {
  small_ccs_df <- data.frame(
    visit_id = c("b", "j", "b"),
    icd9 = c("E8490", "E0000", "E9286")
  )

  # test_ccs_map <- icd9_map_single_ccs[282:284]
  # test_ccs_map <- lapply(test_ccs_map, `[[`, 1)
  #test_ccs_map <- lapply(test_ccs_map, as.character)
  test_ccs_map <- list(`2619` = "E9286",
                       `2620` = "E0000",
                       `2621` = "E8490")
  # this simple map results in the map being the identity matrix

  expected_res <- matrix(byrow = TRUE,
                         data = c(TRUE, FALSE, TRUE,
                                  FALSE, TRUE, FALSE),
                         nrow = 2, ncol = 3,
                         dimnames = list(c("b", "j"),
                                         c("2619", "2620", "2621"))
                         )

  res <- icd9_comorbid_ccs(small_ccs_df, map = test_ccs_map)
  res2 <- icd9_comorbid_ccs(small_ccs_df, map = test_ccs_map, comorbid_fun = icd:::icd9ComorbidShortCpp)
  # compare all three, for development only
  expect_identical(res, expected_res)
  expect_identical(res2, expected_res)
  expect_identical(res, res2)
  # should be one for each (do this way to ignore the empty first group)
  expect_equal(sum(res), length(first_from_each))
})

test_that("ahrq ccs icd 9 is performing correctly", {
  test_df <-
    data.frame(
      visit_id = c("a", "b", "b", "c"),
      icd9 = c("01012", "32341", "83314", "7721"),
      single = c("1", "77", "225", "224"),
      lvl1 = c("1", "6", "16", "15"),
      lvl2 = c("1.1", "6.1", "16.1", "15.7"),
      lvl3 = c("1.1.1", "6.1.2", " ", "15.7.4"),
      lvl4 = c(" ", " ", " ", " "),
      stringsAsFactors = FALSE
    )

  res <- icd9_comorbid_ccs(test_df,  visit_name = "visit_id", icd_name = "icd9")
  # just run this if we find the _alt_ function available, not in production build.
  if (length(grep("icd9Comorbid_alt_MatMul", unclass(lsf.str(envir = asNamespace("icd"), all = TRUE)))) > 0)
    res <- icd9_comorbid_ccs(test_df,  visit_name = "visit_id", icd_name = "icd9", comorbid_fun = icd:::icd9Comorbid_alt_MatMul)
  res2 <- icd9_comorbid_ccs(test_df,  visit_name = "visit_id", icd_name = "icd9", comorbid_fun = icd:::icd9ComorbidShortCpp)

  a_res <- which(sapply(icd9_map_single_ccs, function(y) "01012" %in% y))
  b_res1 <- which(sapply(icd9_map_single_ccs, function(y) "32341" %in% y))
  b_res2 <- which(sapply(icd9_map_single_ccs, function(y) "83314" %in% y))
  c_res <- which(sapply(icd9_map_single_ccs, function(y) "7721" %in% y))

  # build an unnamed matrix with the right flags set
  manual_res <- matrix(FALSE, nrow = 3, ncol = 284)
  manual_res[1, a_res] <- TRUE
  manual_res[2, b_res1] <- TRUE
  manual_res[2, b_res2] <- TRUE
  manual_res[3, c_res] <- TRUE
  expect_equivalent(manual_res, res2)
  expect_equivalent(manual_res, res)
  expect_identical(res, res2)
  expect_true(all(mapply(function(x, y) res[x, y], test_df$visit_id, test_df$single)))
  expect_equal(dim(res), c(3, 284))

  expect_error(icd9_comorbid_ccs(test_df, visit_name = "visit_id", icd_name = "icd9", single = FALSE))
  expect_error(icd9_comorbid_ccs(test_df, visit_name = "visit_id", icd_name = "icd9", single = FALSE, lvl = "a"))

  res <- icd9_comorbid_ccs(test_df, visit_name = "visit_id", icd_name = "icd9", single = FALSE, lvl = 1)
  expect_true(all(mapply(function(x, y) res[x, y], test_df$visit_id, test_df$lvl1)))
  expect_equal(dim(res), c(3, 18))

  res <- icd9_comorbid_ccs(test_df, visit_name = "visit_id", icd_name = "icd9", single = FALSE, lvl = 2)
  expect_true(all(mapply(function(x, y) res[x, y], test_df$visit_id, test_df$lvl2)))
  expect_equal(dim(res), c(3, 136))

  res <- icd9_comorbid_ccs(test_df, visit_name = "visit_id", icd_name = "icd9",  single = FALSE,  lvl = 3)
  expect_true(all(mapply(function(x, y) res[x, y], test_df$visit_id, test_df$lvl3)))
  expect_equal(dim(res), c(3, 367))

  res <- icd9_comorbid_ccs(test_df, visit_name = "visit_id", icd_name = "icd9", single = FALSE, lvl = 4)
  expect_true(all(mapply(function(x, y) res[x, y], test_df$visit_id, test_df$lvl4)))
  expect_equal(dim(res), c(3, 209))
})
