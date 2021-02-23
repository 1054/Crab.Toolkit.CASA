# 
# Run this code inside CASA as:
#     execfile('')
# or
#     exec(open('').read())
# 

import os, sys, re, glob, shutil
sys.path.append(os.path.expanduser('~/Cloud/Github/Crab.Toolkit.CASA/lib/python'))
from dzliu_clean_utils import get_datacolumn
from taskinit import casalog, tb


# User-defined parameters
DataSet_dirs = sorted(glob.glob('Level_2_Calib/DataSet_*/calibrated/calibrated.ms'))
DataSet_names = [os.path.basename(os.path.dirname(os.path.dirname((t)))) for t in DataSet_dirs]
Target_name = 'GOODS-S'
Split_width = 2 # 31.25MHz
Split_spw = '5,7,9,11' # 
output_concat_vis = 'Level_3_Split/DataSet_Merged/split_GOODS-S_sci.ms'


# def print2
def print2(message):
    print(message)
    casalog.post(message, 'INFO')


# Start processing
print2('getcwd: %s'%(os.getcwd()))

list_of_vis_to_concat = []

for DataSet_dir, DataSet_name in list(zip(DataSet_dirs, DataSet_names)):

    vis = os.path.abspath(DataSet_dir)
    field = Target_name
    outputvis = 'split_'+Target_name+'_sci.ms'
    width = Split_width
    spw = Split_spw
    current_dir = os.getcwd()
    list_of_vis_to_concat.append('Level_3_Split/'+DataSet_name+'/'+outputvis)
    
    # prepare output dir
    if not os.path.isdir('Level_3_Split/'+DataSet_name):
        os.makedirs('Level_3_Split/'+DataSet_name)
    os.chdir('Level_3_Split/'+DataSet_name)
    print2('chdir: %s'%(os.getcwd()))
    
    # run listobs
    if not os.path.isfile(vis+'.listobs.txt'):
        print2('Running listobs for vis %r and writing into %r'%(vis, vis+'.listobs.txt'))
        listobs(vis=vis, listfile=vis+'.listobs.txt')
        #raise NotImplementedError()

    # run split
    if not os.path.isdir(outputvis):
        split_params = {}
        split_params['vis'] = vis
        split_params['outputvis'] = outputvis
        split_params['field'] = field
        split_params['timebin'] = '30s'
        split_params['keepflags'] = True
        split_params['keepmms'] = False
        split_params['datacolumn'] = get_datacolumn(vis)
        split_params['width'] = width
        split_params['spw'] = spw
        print2('Splitting ms ...')
        print2('split('+', '.join("{!s}={!r}".format(k, split_params[k]) for k in split_params.keys())+')')
        split(**split_params)
        shutil.copy2('split.last', outputvis+'.split.last')
        if not os.path.isdir(outputvis):
            raise Exception('Error! Failed to split the ms.')

    # cd back
    os.chdir(current_dir)    
    print2('chdir: %s'%(os.getcwd()))    



# concat
if not os.path.isdir(output_concat_vis):
    concat_params = {}
    concat_params['vis'] = list_of_vis_to_concat
    concat_params['concatvis'] = output_concat_vis
    concat_params['freqtol'] = '20MHz'
    concat_params['copypointing'] = True #<TODO># 
    print2('Concatenating ms ...')
    print2('concat('+', '.join("{!s}={!r}".format(k, concat_params[k]) for k in concat_params.keys())+')')
    concat(**concat_params)
    shutil.copy2('concat.last', output_concat_vis+'.concat.last')
    if not os.path.isdir(output_concat_vis):
        raise Exception('Error! Failed to concatenate the ms.')










