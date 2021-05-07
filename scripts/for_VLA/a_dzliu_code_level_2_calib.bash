#!/bin/bash
# 
set -e

source ~/Cloud/Github/AlmaCosmos/Software/SETUP.bash
source ~/Softwares/CASA/SETUP.bash 5.0.0
#source ~/Softwares/GILDAS/SETUP.bash
#source ~/Cloud/Github/Crab.Toolkit.PdBI/SETUP.bash


# First fix scan xml
#./a_dzliu_code_level_2_calib_task_1_fix_scan_xml.py


# Then make data directory structure
#alma_archive_make_data_dirs_with_meta_table.py meta_data_table.txt


# Then cd into each directory and run pipeline calibration
DataSet_list=($(ls -1d Level_2_Calib/DataSet_* | sort -V))
for (( i = 0; i < ${#DataSet_list[@]}; i++ )); do
    DataSet_dir=$(basename "${DataSet_list[i]}")
    cd Level_2_Calib/$DataSet_dir/calibrated/
    if [[ ! -d calibrated.ms ]] && [[ ! -L calibrated.ms ]]; then
        # 
        ./scriptForDatasetRecalibration.sh
        # 
        MeasurementSet_list=($(find . -maxdepth 1 -type d -name "*.ms" -a -not -name "calibrators.ms" -a -not -name "split*.ms"))
        if [[ ${#MeasurementSet_list[@]} -gt 1 ]]; then
            echo "Error! There are more than one(${#MeasurementSet_list[@]}) ms under \"$(pwd)\"!"
            echo "${MeasurementSet_list[@]}"
            exit 255
        fi
        if [[ ${#MeasurementSet_list[@]} -ge 1 ]]; then
            echo "ln -fsT ${MeasurementSet_list[0]} calibrated.ms"
            ln -fsT "${MeasurementSet_list[0]}" "calibrated.ms"
        fi
    fi
    cd ../../../
done


# Then cd into each directory again and link measurement sets to a common named link "calibrated.ms"
#./a_dzliu_code_level_2_calib_task_2_make_links.bash



