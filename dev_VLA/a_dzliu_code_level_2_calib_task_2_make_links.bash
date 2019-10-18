#!/bin/bash
# 


# cd into each directory and check calibrated ms and make links
DataSet_list=($(ls -1d Level_2_Calib/DataSet_* | sort -V))
for (( i = 0; i < ${#DataSet_list[@]}; i++ )); do
    DataSet_dir=$(basename "${DataSet_list[i]}")
    cd Level_2_Calib/$DataSet_dir/calibrated/
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
    cd ../../../
done


