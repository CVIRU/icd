#!/bin/bash
set -uo pipefail
IFS=$'\n\t'

ICD_HOME=${ICD_HOME:-"$HOME/rprojects/icd"}
INSTR_ATSTART="no"
VALGRIND_CMD="valgrind --tool=callgrind --simulate-cache=yes --instr-atstart=$INSTR_ATSTART --separate-threads=no"

crashscript="${ICD_HOME:-/tools/crash-scripts/eigen.R}"

function valgrind-cmd () {
	pushd "$ICD_HOME"
	R --vanilla --slave -d "$VALGRIND_CMD" -e \'$RCODE\'
	popd
}

function valgrind-script () {
	pushd "$ICD_HOME"
	R --vanilla --slave -d "$VALGRIND_CMD" < "${crashscript}"
	popd
}

