#!/bin/bash
# 

source ~/Software/CASA/SETUP.bash 5.0.0

script_name=$(basename "${BASH_SOURCE[0]}" | perl -p -e 's/\.bash$//g')
casa --nogui --log2term -c "exec(compile(open(\"${script_name}.py\").read(),\"${script_name}.py\",\"exec\"))"

