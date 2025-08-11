#!/bin/bash
# 

#source ~/Software/CASA/SETUP.bash 6.5.5 # uvcontsub changed! no fitspw but fitspec, no excludechans, use uvcontsub_old
source ~/Software/CASA/SETUP.bash 5.7.2

script_name=$(basename ${BASH_SOURCE[0]/.bash//})

echo "casa --nogui --log2term -c \"execfile('${script_name}.py')\" 2>/dev/null | tee \"log_${script_name}.txt\""
casa --nogui --log2term -c "execfile('${script_name}.py')" 2>/dev/null | tee "log_${script_name}.txt"

