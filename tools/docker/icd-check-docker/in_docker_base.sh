#!/bin/bash
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

set -euo pipefail
IFS=$'\n\t'

# echo expanded commands
set -x

#https://stackoverflow.com/questions/14719349/error-c-stack-usage-is-too-close-to-the-limit#14719448
# C stack limit problems hopefully fixed by this. Can also be set to "unlimited" Default is 8192k on my pc.
# do this inside and outside the docker container

old_ulimit=$(ulimit -s)
ulimit -s unlimited

function finish {
        ulimit -s "$old_ulimit"
}
trap finish EXIT

pushd /tmp

if [ -f /etc/os-release ]; then
    source /etc/os-release
    OS=$NAME
elif type lsb_release >/dev/null 2>&1; then
    OS=$(lsb_release -si)
elif [ -f /etc/lsb-release ]; then
    source /etc/lsb-release
    OS=$DISTRIB_ID
elif [ -f /etc/debian_version ]; then
    OS=Debian
else
    OS=$(uname -s)
fi

if [[ "$OS" == "Ubuntu" || "$OS" == "Debian" ]]; then
  echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
  apt-get update -qq || true
  apt-get dist-upgrade -qq -y || true
  apt-get install -y -qq git libssl-dev libxml2-dev curl libcurl4-openssl-dev unixodbc-dev qpdf pandoc pandoc-citeproc # libssh2-1-dev (optional for git but has debian version problem at least in April 2017)
elif [[ "$OS" == "Fedora" ]]; then
  dnf install -y git
fi

echo "R CMD: using '${R_CMD:-RD}'"
echo "R CMD ERR: using '${R_CMD_ERR:-RD}'"

R_CMD_ERR="${R_CMD}script -e 'cat()'"

# which libasan library?
if [ -e /usr/lib/llvm-4.0/lib/clang/4.0.1/lib/linux/libclang_rt.asan-x86_64.so ]; then
  if LD_PRELOAD="/usr/lib/llvm-4.0/lib/clang/4.0.1/lib/linux/libclang_rt.asan-x86_64.so" ${R_CMD_ERR}; then
    export LD_PRELOAD="/usr/lib/llvm-4.0/lib/clang/4.0.1/lib/linux/libclang_rt.asan-x86_64.so"
  fi
elif [ -e /usr/lib/x86_64-linux-gnu/libasan.so.4 ]; then
  if LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libasan.so.4" ${R_CMD_ERR}; then
    export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libasan.so.4"
  fi
elif [ -e /usr/lib/x86_64-linux-gnu/libasan.so.3 ]; then
  if LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libasan.so.3" ${R_CMD_ERR}; then
    export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libasan.so.3"
  fi
elif [ -e /usr/lib/x86_64-linux-gnu/libasan.so.2 ]; then
  if LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libasan.so.2" ${R_CMD_ERR}; then
    export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libasan.so.2"
  fi
else
  echo "Cannot find libasan in /usr/lib/x86_64-linux-gnu"
fi

# tolerate R_CMD unset or empty, and default to RD if empty or unset:
# TODO: actually, RD is not available in all docker images, e.g. most basic rocker/tidyverse, verse, etc.

echo "checking ${R_CMD} exists"
if ! command -v ${R_CMD:-RD} &>/dev/null; then
  echo "setting R_CMD to R"
  R_CMD=R
fi

# these are always checked for, so we don't care which R is installed. We also need to re-install some packages, for some reason unclear to me: https://github.com/rocker-org/r-devel-san-clang/issues/12
for pkg in testthat checkmate RODBC xml2 Rcpp stringi knitr rmarkdown microbenchmark; do
ASAN_OPTIONS=abort_on_error=0,detect_leaks=0 ${R_CMD}script -e "install.packages(\"${pkg}\")"
done


# We need to build based on the directory name, not the package name. Don't build vignettes now: they get built during check anyway.
# $R_CMD CMD build --no-build-vignettes $GITHUB_REPO 
popd # out of /tmp

