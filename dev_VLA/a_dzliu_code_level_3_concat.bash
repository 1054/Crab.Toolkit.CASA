#!/bin/bash
# 

source ~/Softwares/CASA/SETUP.bash 5.0.0
source ~/Softwares/GILDAS/SETUP.bash
#source ~/Cloud/Github/Crab.Toolkit.PdBI/SETUP.bash


# Prepare output directory
if [[ ! -d "Level_3_Concat" ]]; then
    mkdir "Level_3_Concat"
    if [[ ! -d "Level_3_Concat" ]]; then
        echo "Error! Failed to make \"Level_3_Concat\" directory!"
        exit 255
    fi
fi

# 
# Run script to match source names and concat
./a_dzliu_code_level_3_concat_task_1_casa_concat.py


