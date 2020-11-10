#!/usr/bin/env bash
# 


#Crab_BIN_SETUP_SCRIPT=$(dirname "${BASH_SOURCE[0]}")/bin/bin_setup.bash

#source "$Crab_BIN_SETUP_SCRIPT" -check casa-ms-split casa-ms-concat pdbi-uvt-go-splitpolar pdbi-uvt-go-import-uvfits pdbi-uvt-go-average pdbi-uvt-go-uvfit

export Crab_Toolkit_CASA_PATH=$(dirname $(perl -MCwd -e 'print Cwd::abs_path shift' "${BASH_SOURCE[0]}"))


