Resubmitting after removing accidental unicode character inclusion. Thank you.

# Test environments
* Ubuntu 17.10 R 3.4.3, gcc 7.2, clang 5, clang-trunk
* Ubuntu 14.04 (on travis-ci) R-devel, gcc
* Debian in docker with R-devel, ASAN+UBSAN GCC 7.2, clang 5 & trunk
* Appveyor Windows Server 2012 R2 x64, R 3.4.3 32 and 64 bit

# R CMD check results

There is one note:

checking installed package size ... NOTE
  installed size is 10.3Mb
  sub-directories of 1Mb or more:
    R      1.5Mb
    data   2.1Mb
    doc    3.2Mb
    libs   3.0Mb
