# -*- coding: utf-8 -*-
# 
"""Utilities for analyzing u-v visibility data within CASA.

Notes
-----
Functions in this code must be run in CASA.  

Functions
---------
This code contains following functions (not a complete list):

- TBD

Last updates
------------
- 2021-08-16 start

Example
-------
Example commands to run this code (in CASA environment)::

    import os, sys, glob
    sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
    from dzliu_uvdata_utils import plot_uvdist_vs_amp
    plot_uvdist_vs_amp(vis, plotfile)

"""

from __future__ import print_function
import os, sys, re, json, copy, time, datetime, shutil
import numpy as np


# 
# class json_info_dict_encoder(json.JSONEncoder)
# -- for json dump and load 
# -- https://stackoverflow.com/questions/27050108/convert-numpy-type-to-python
# 
class json_info_dict_encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(json_info_dict_encoder, self).default(obj)


# 
# Function to plot uvdist vs amp
# 
def plot_uvdist_vs_amp(vis, plotfile):
    ..
    
    





