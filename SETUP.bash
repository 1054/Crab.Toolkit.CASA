#!/usr/bin/env bash
# 


#Crab_BIN_SETUP_SCRIPT=$(dirname "${BASH_SOURCE[0]}")/bin/bin_setup.bash

#source "$Crab_BIN_SETUP_SCRIPT" -check casa-ms-split casa-ms-concat pdbi-uvt-go-splitpolar pdbi-uvt-go-import-uvfits pdbi-uvt-go-average pdbi-uvt-go-uvfit

export Crab_Toolkit_CASA_PATH=$(dirname $(perl -MCwd -e 'print Cwd::abs_path shift' "${BASH_SOURCE[0]}"))


# write startup file
# see https://casaguides.nrao.edu/index.php/Analysis_Utilities
if [[ ! -d "$HOME/.casa" ]]; then
    mkdir -p "$HOME/.casa"
fi

do_write_file=1
if [[ -f "$HOME/.casa/init.py" ]]; then
    if [[ $(grep "${Crab_Toolkit_CASA_PATH}/lib/python" "$HOME/.casa/init.py" | wc -l) -ge 1 ]]; then
        do_write_file=0
    fi
fi
if [[ $do_write_file -eq 1 ]]; then
    echo "" >> "$HOME/.casa/init.py"
    echo "import sys" >> "$HOME/.casa/init.py"
    echo "sys.path.append(\"${Crab_Toolkit_CASA_PATH}/lib/python\")" >> "$HOME/.casa/init.py"
    echo "" >> "$HOME/.casa/init.py"
fi

do_write_file=1
if [[ -f "$HOME/.casa/startup.py" ]]; then
    if [[ $(grep "${Crab_Toolkit_CASA_PATH}/lib/python" "$HOME/.casa/startup.py" | wc -l) -ge 1 ]]; then
        do_write_file=0
    fi
fi
if [[ $do_write_file -eq 1 ]]; then
    echo "" >> "$HOME/.casa/startup.py"
    echo "import sys" >> "$HOME/.casa/startup.py"
    echo "sys.path.append(\"${Crab_Toolkit_CASA_PATH}/lib/python\")" >> "$HOME/.casa/startup.py"
    echo "" >> "$HOME/.casa/startup.py"
fi
        

