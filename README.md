# Crab.Toolkit.CASA

This repository contains some python libraries for ALMA data processing in CASA software.

Below are some examples. We assume that this repo has been downloaded into a local folder `~/Cloud/Github/Crab.Toolkit.CASA`. 

## For Imaging

We can run following commands in CASA for easy imaging of the ALMA data. These commands can also be written in a script which can then be executed with the `execfile` command in CASA. 
```
# Run this code in CASA
import os, sys 
sys.path.insert(1, os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
import dzliu_clean
vis = '/path/to/your/split/science/target/data.ms'
output_image = '/path/to/the/output/image' # without suffix
galaxy_name = 'your_target_name'
dzliu_clean.dzliu_clean(vis, output_image = output_image, galaxy_name = galaxy_name, galaxy_redshift = None, 
                        make_line_cube = True, line_name = 'cube', line_velocity = -1, line_velocity_width = -1, 
                        line_velocity_resolution = 20.0, line_clean_threshold = 2.0, max_imsize = 1500)
                        # here we set the line_velocity_resolution to 20.0 km/s, but the actual output channel width
                        # will be an integer factor of the original channel width as close to 20.0 km/s as possible. 
```

Or we can do more things with the utils library. Note that the input visibility data measurement set (ms) should be a split ms with only one science target (field) and one spectral window (spw).  
```
# Run this code in CASA
import os, sys 
sys.path.insert(1, os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
from dzliu_clean_utils import (get_datacolumn, get_mosaic_imsize_and_phasecenter, get_synbeam_and_imcell, 
                               get_spw_for_spectral_line, 
                               cleanup_tclean_products, apply_pbcor_to_tclean_image, export_tclean_products_as_fits_files)
vis = '/path/to/your/split/science/target/data.ms'
datacolumn = get_datacolumn(vis)
beam, cell = get_synbeam_and_imcell(vis)
spw = get_spw_for_spectral_line(vis, redshift=None, rest_freq_GHz=None, line_width_kms=None)
......
```



## For combining uvfits from different projects






