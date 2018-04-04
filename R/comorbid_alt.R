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

#' Comorbidities from ICD-10 codes using matrix multiplication
#'
#' Use the fact that we only need to look up comorbidities for those codes which
#' are in the input data frame, and map. Difference with ICD-9 codes is that we
#' only have 'parent' codes in the map.
#'
#' 1. get list of unique ICD-10 codes in the source data
#'
#' 2. make a new comorbidity map containing those codes where those codes'
#' parents appear in the given map
#'
#' 3. use the new common comorbidity code using openmp eigen matrix multiply
#' with that data.
#' @keyword internal
icd10_comorbid_matmul <- function(
  x,
  map,
  visit_name = NULL,
  icd_name = get_icd_name(x),
  short_code = icd_guess_short(x, icd_name = icd_name),
  short_map = icd_guess_short(map),
  return_df = FALSE, ...) {

  uniq_codes <- unique(x[[icd_name]])


}

icd10_comorbid_parent_search_str <- function(
  x,
  map,
  visit_name = NULL,
  icd_name = get_icd_name(x),
  short_code = icd_guess_short(x, icd_name = icd_name),
  short_map = icd_guess_short(map),
  return_df = FALSE, ...) {

  if (!short_code)
    x[[icd_name]] <- icd_decimal_to_short.icd10(x[[icd_name]])

  icd_codes <- x[[icd_name]]

  # for each icd code
  just_cmb <- vapply(icd_codes, FUN.VALUE = logical(length(map)), FUN = function(y) {
    # look it up in each comorbidity, but SOMEDAY: once we have a comorbidity for
    # one patient, we don't need to search within it again

    char_count <- nchar(as.character(y)):3
    vapply(names(map), FUN.VALUE = logical(1),
           FUN = function(cmb) {
             # and if not found, slice off last char of test string
             for (n in char_count) {
               if (!is.na(match(substr(y, 1, n), map[[cmb]])))
                 return(TRUE)
             }
             FALSE
           })
  })

  res <- aggregate(x = t(just_cmb), by = x[visit_name], FUN = any)
  if (return_df)
    return(res)

  out <- as.matrix(res[-1])
  rownames(out) <- res[[1]]
  out
}

icd10_comorbid_parent_search_all <- function(x,
                                             map,
                                             visit_name = NULL,
                                             icd_name = get_icd_name(x),
                                             short_code = icd_guess_short(x, icd_name = icd_name),
                                             short_map = icd_guess_short(map),
                                             return_df = FALSE, ...) {

  if (!short_code)
    x[[icd_name]] <- icd_decimal_to_short.icd10(x[[icd_name]])

  icd_codes <- x[[icd_name]]

  # for each icd code
  just_cmb <- vapply(icd_codes, FUN.VALUE = logical(length(map)), FUN = function(y) {
    # look it up in each comorbidity, but SOMEDAY: once we have a comorbidity for
    # one patient, we don't need to search within it again


    char_count <- nchar(as.character(y)):3
    vapply(names(map), FUN.VALUE = logical(1),
           FUN = function(cmb) {
             # and if not found, slice off last char of test string
             perms_to_match <- vapply(char_count:3, substr, x = y, start = 1, FUN.VALUE = character(1))
             any(match(perms_to_match, map[[cmb]], nomatch = 0L) != 0)
           })
  })

  res <- aggregate(x = t(just_cmb), by = x[visit_name], FUN = any)
  if (return_df)
    return(res)

  out <- as.matrix(res[-1])
  rownames(out) <- res[[1]]
  out
}

icd10_comorbid_parent_search_no_loop <- function(x,
                                                 map,
                                                 visit_name = NULL,
                                                 icd_name = get_icd_name(x),
                                                 short_code = icd_guess_short(x, icd_name = icd_name),
                                                 short_map = icd_guess_short(map),
                                                 return_df = FALSE, ...) {

  if (!short_code)
    x[[icd_name]] <- icd_decimal_to_short.icd10(x[[icd_name]])

  icd_codes <- x[[icd_name]]

  # for each icd code
  just_cmb <- vapply(icd_codes, FUN.VALUE = logical(length(map)), FUN = function(y) {
    # look it up in each comorbidity, but SOMEDAY: once we have a comorbidity for
    # one patient, we don't need to search within it again

    vapply(names(map), FUN.VALUE = logical(1),
           FUN = function(cmb) {
             # instead of loop, just declare the substring length
             if (!is.na(match(substr(y, 1, 10), map[[cmb]]))) return(TRUE)
             if (!is.na(match(substr(y, 1, 9), map[[cmb]]))) return(TRUE)
             if (!is.na(match(substr(y, 1, 8), map[[cmb]]))) return(TRUE)
             if (!is.na(match(substr(y, 1, 7), map[[cmb]]))) return(TRUE)
             if (!is.na(match(substr(y, 1, 6), map[[cmb]]))) return(TRUE)
             if (!is.na(match(substr(y, 1, 5), map[[cmb]]))) return(TRUE)
             if (!is.na(match(substr(y, 1, 4), map[[cmb]]))) return(TRUE)
             if (!is.na(match(substr(y, 1, 3), map[[cmb]]))) return(TRUE)
             FALSE
           })
  })

  res <- aggregate(x = t(just_cmb), by = x[visit_name], FUN = any)
  if (return_df)
    return(res)

  out <- as.matrix(res[-1])
  rownames(out) <- res[[1]]
  out
}

icd10_comorbid_parent_search_orig <- function(x,
                                              map,
                                              visit_name = NULL,
                                              icd_name = get_icd_name(x),
                                              short_code = icd_guess_short(x, icd_name = icd_name),
                                              short_map = icd_guess_short(map),
                                              return_df = FALSE, ...) {

  if (!short_code)
    x[[icd_name]] <- icd_decimal_to_short.icd10(x[[icd_name]])

  icd_codes <- x[[icd_name]]

  # for each icd code
  just_cmb <- vapply(icd_codes, FUN.VALUE = logical(length(map)), FUN = function(y) {
    # look it up in each comorbidity, but SOMEDAY: once we have a comorbidity for
    # one patient, we don't need to search within it again


    char_count <- nchar(as.character(y)):3
    vapply(names(map), FUN.VALUE = logical(1),
           FUN = function(cmb) {
             # and if not found, slice off last char of test string
             for (n in char_count) {
               if (!is.na(match(substr(y, 1, n), map[[cmb]])))
                 return(TRUE)
             }
             FALSE
           })
  })

  res <- aggregate(x = t(just_cmb), by = x[visit_name], FUN = any)
  if (return_df)
    return(res)

  out <- as.matrix(res[-1])
  rownames(out) <- res[[1]]
  out
}

icd10_comorbid_no_parent_search <- function(x,
                                            map,
                                            visit_name = NULL,
                                            icd_name = get_icd_name(x),
                                            short_code = icd_guess_short(x, icd_name = icd_name),
                                            short_map = icd_guess_short(map[[1]]),
                                            return_df = FALSE, ...) {

  # confirm class is ICD-9 so we dispatch correctly. The class may not be set if
  # the S3 method was called directly.
  if (!is.icd10(x[[icd_name]])) x[[icd_name]] <- icd10(x[[icd_name]]) # no as.icd10
  icd_comorbid_common(x, map, visit_name, icd_name,
                      short_code, short_map, return_df, ...)

}
