R --slave --vanilla -d lldb -e "library(icd); library(magrittr); icd_comorbid_ahrq(vermont_dx %>% icd_wide_to_long, comorbid_fun = icd:::icd9Comorbid_alt_MatMul)"

