
# RUN THIS SCRIPT INSIDE CASA

# Read *.ms or *.uvfits data file from the system variables.
# Run CASA tclean and output image cube and other files.

import os, sys
import numpy as np
CASA_COMMANDLINE_VIS = os.getenv('CASA_COMMANDLINE_VIS')
CASA_COMMANDLINE_OUT = os.getenv('CASA_COMMANDLINE_OUT')
CASA_COMMANDLINE_FIELD = os.getenv('CASA_COMMANDLINE_FIELD')
CASA_COMMANDLINE_SYSPATH = os.getenv('CASA_COMMANDLINE_SYSPATH')
CASA_COMMANDLINE_OVERWRITE = os.getenv('CASA_COMMANDLINE_OVERWRITE', 'False').lower() in ('true', '1', 't')

if CASA_COMMANDLINE_SYSPATH not in sys.path:
    sys.path.append(CASA_COMMANDLINE_SYSPATH)
    sys.path.append(os.path.join(CASA_COMMANDLINE_SYSPATH, 'analysis_scripts'))

#from dzliu_clean_utils import (get_datacolumn, get_synbeam_and_imcell, get_mstransform_params_for_spectral_line, 
#                               apply_pbcor_to_tclean_image, export_tclean_products_as_fits_files, cleanup_tclean_products)

import dzliu_clean

dzliu_clean.dzliu_clean(
    CASA_COMMANDLINE_VIS, 
    output_image = CASA_COMMANDLINE_OUT, 
    make_line_cube = False, 
    make_continuum = True, 
    overwrite = CASA_COMMANDLINE_OVERWRITE, 
    )













