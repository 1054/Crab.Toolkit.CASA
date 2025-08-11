#!/bin/bash
# 

find . -maxdepth 1 -type d -name "input_data*" -exec rm -rf {} \;

find run_tclean_perplanebeam -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} \;

find run_tclean_perplanebeam_chw15kms -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} \;

[ -d combined.ms ] && rm -rf combined.ms*






