#!/bin/bash
# 
set -e

RA="00:00:00"
Dec="00:00:00"


# Make working dir
if [[ ! -d "Level_4_Clean/run_tclean_pointing6/run_gildas_uvfit" ]]; then
mkdir -p "Level_4_Clean/run_tclean_pointing6/run_gildas_uvfit"
fi

cd "Level_4_Clean/run_tclean_pointing6/"



# First detect bright continuum sources with PyBDSF
if [[ ! -f Output_Blind_Extraction_Photometry_PyBDSM/pointing6_cont_clean.image/pybdsm_cat.fits ]]; then
    ~/Cloud/Github/AlmaCosmos/Software/AlmaCosmos_Photometry_Blind_Extraction_PyBDSM_mod.py pointing6_cont_clean.image.fits -thresh_rms 6.0
fi
if [[ ! -f Output_Blind_Extraction_Photometry_PyBDSM/pointing6_cont_clean.image/pybdsm_cat.ascii.txt ]]; then
    topcat -stilts tpipe in=Output_Blind_Extraction_Photometry_PyBDSM/pointing6_cont_clean.image/pybdsm_cat.fits \
                        cmd='keepcols "RA DEC"' ofmt=ascii \
                        out=Output_Blind_Extraction_Photometry_PyBDSM/pointing6_cont_clean.image/pybdsm_cat.ascii.txt
fi

list_bright_sources_RA=($(cat Output_Blind_Extraction_Photometry_PyBDSM/pointing6_cont_clean.image/pybdsm_cat.ascii.txt | grep -v '^#' | sed -e 's/^ *//g' | cut -d ' ' -f 1))
list_bright_sources_Dec=($(cat Output_Blind_Extraction_Photometry_PyBDSM/pointing6_cont_clean.image/pybdsm_cat.ascii.txt | grep -v '^#' | sed -e 's/^ *//g' | cut -d ' ' -f 2))


# cd working dir
cd "run_gildas_uvfit"


# Prepare to run GILDAS uvfit
#source ~/Cloud/Github/Crab.Toolkit.PdBI/SETUP.bash
#source ~/Softwares/GILDAS/SETUP.bash
#source ~/Softwares/CASA/SETUP.bash 5.0.0


# convert ms to uvt
script_file="a_dzliu_code_step_1_convert_ms_to_uvt.bash"
if [[ ! -f split_pointing6_spw0_width1_SP.uvt ]]; then
echo "#!/bin/bash" > "$script_file"
echo "#" >> "$script_file"
echo "source \$HOME/Cloud/Github/Crab.Toolkit.PdBI/SETUP.bash" >> "$script_file"
echo "source \$HOME/Softwares/CASA/SETUP.bash 5.0.0" >> "$script_file"
echo "source \$HOME/Softwares/GILDAS/SETUP.bash" >> "$script_file"
echo "" >> "$script_file"
echo "casa-ms-split -vis ../pointing6_HI21cm.ms -step split exportuvfits gildas -width 1" >> "$script_file"
echo "" >> "$script_file"
chmod +x "$script_file"
cat "$script_file"
./"$script_file"
fi


script_file="a_dzliu_code_step_2_merge_spws.bash"
if [[ ! -f merged.uvt ]]; then
echo "#!/bin/bash" > "$script_file"
echo "#" >> "$script_file"
echo "source \$HOME/Cloud/Github/Crab.Toolkit.PdBI/SETUP.bash" >> "$script_file"
echo "source \$HOME/Softwares/CASA/SETUP.bash 5.0.0" >> "$script_file"
echo "source \$HOME/Softwares/GILDAS/SETUP.bash" >> "$script_file"
echo "" >> "$script_file"
echo "if [[ \$(ls split_pointing6_spw*_width1_SP.uvt | wc -l) -gt 1 ]]; then" >> "$script_file"
echo "pdbi-uvt-go-merge -name split_pointing6_spw*_width1_SP.uvt -out merged.uvt" >> "$script_file"
echo "else" >> "$script_file"
echo "cp split_pointing6_spw0_width1_SP.uvt merged.uvt" >> "$script_file"
echo "fi" >> "$script_file"
echo "" >> "$script_file"
echo "echo \"header merged.uvt\" | mapping -nw -nl > merged.uvt.header.txt" >> "$script_file"
echo "" >> "$script_file"
chmod +x "$script_file"
cat "$script_file"
./"$script_file"
fi


# Prepare to subtract bright continuum sources
script_file="a_dzliu_code_step_3_subtract_bright_continuum_sources.bash"
if [[ ! -f merged_go_uvfit_residual.uvt ]]; then
echo "#!/bin/bash" > "$script_file"
echo "#" >> "$script_file"
echo "source \$HOME/Cloud/Github/Crab.Toolkit.PdBI/SETUP.bash" >> "$script_file"
echo "source \$HOME/Softwares/CASA/SETUP.bash 5.0.0" >> "$script_file"
echo "source \$HOME/Softwares/GILDAS/SETUP.bash" >> "$script_file"
echo "" >> "$script_file"
echo "pdbi-uvt-go-uvfit -name merged \\" >> "$script_file"
for (( i = 0; i < ${#list_bright_sources_RA[@]}; i++ )); do
echo "    -radec ${list_bright_sources_RA[i]} ${list_bright_sources_Dec[i]} -fixedpos -subtract \\" >> "$script_file"
done
echo "    -out merged_go_uvfit \\" >> "$script_file"
echo "    -res merged_go_uvfit_residual" >> "$script_file"
echo "" >> "$script_file"
chmod +x "$script_file"
cat "$script_file"
./"$script_file"
fi



# Shift phasecenter to galaxy center
script_file="a_dzliu_code_step_4_shift_phase_center_to_galaxy_center.bash"
if [[ ! -f uvshifted.uvt ]]; then
echo "#!/bin/bash" > "$script_file"
echo "#" >> "$script_file"
echo "source \$HOME/Cloud/Github/Crab.Toolkit.PdBI/SETUP.bash" >> "$script_file"
echo "source \$HOME/Softwares/CASA/SETUP.bash 5.0.0" >> "$script_file"
echo "source \$HOME/Softwares/GILDAS/SETUP.bash" >> "$script_file"
echo "" >> "$script_file"
echo "pdbi-uvt-go-shift -name merged_go_uvfit_residual.uvt -radec ${RA} ${Dec} -out uvshifted.uvt" >> "$script_file"
echo "" >> "$script_file"
chmod +x "$script_file"
cat "$script_file"
./"$script_file"
fi



# Extract fixed-position spectrum
script_file="a_dzliu_code_step_5_extract_fixed_position_point_source_spectrum.bash"
if [[ ! -f uvshifted_go_uvfit_target.result.obj_1.spectrum.pdf ]]; then
echo "#!/bin/bash" > "$script_file"
echo "#" >> "$script_file"
echo "source \$HOME/Cloud/Github/Crab.Toolkit.PdBI/SETUP.bash" >> "$script_file"
echo "source \$HOME/Softwares/CASA/SETUP.bash 5.0.0" >> "$script_file"
echo "source \$HOME/Softwares/GILDAS/SETUP.bash" >> "$script_file"
echo "" >> "$script_file"
echo "pdbi-uvt-go-uvfit -name uvshifted -offset 0 0 -fixedpos -out uvshifted_go_uvfit_target" >> "$script_file"
echo "" >> "$script_file"
chmod +x "$script_file"
cat "$script_file"
./"$script_file"
fi








