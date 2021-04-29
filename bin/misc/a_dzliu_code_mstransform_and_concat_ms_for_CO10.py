
# Run this script in CASA

import os, sys, re, shutil
from collections import OrderedDict

list_of_ms = [\
'../VLA_10A-157_AB1347_PI_Baker/Level_2_Calib/DataSet_01/calibrated/calibrated.ms',
'../VLA_10A-157_AB1347_PI_Baker/Level_2_Calib/DataSet_02/calibrated/calibrated.ms',
'../VLA_10A-157_AB1347_PI_Baker/Level_2_Calib/DataSet_03/calibrated/calibrated.ms',
'../VLA_10A-157_AB1347_PI_Baker/Level_2_Calib/DataSet_04/calibrated/calibrated.ms',
'../VLA_10A-157_AB1347_PI_Baker/Level_2_Calib/DataSet_05/calibrated/calibrated.ms',
'../VLA_10C-223_PI_Sharon/Level_2_Calib/DataSet_01/calibrated/calibrated.ms', 
'../VLA_12A-399_PI_Sharon/Level_2_Calib/DataSet_01/calibrated/calibrated.ms', 
'../VLA_12A-399_PI_Sharon/Level_2_Calib/DataSet_02/calibrated/calibrated.ms', 
'../VLA_12A-399_PI_Sharon/Level_2_Calib/DataSet_03/calibrated/calibrated.ms', 
'../VLA_12A-399_PI_Sharon/Level_2_Calib/DataSet_04/calibrated/calibrated.ms', 
'../VLA_12A-399_PI_Sharon/Level_2_Calib/DataSet_05/calibrated/calibrated.ms', 
'../VLA_14A-371_PI_Carilli/Level_2_Calib/DataSet_01/calibrated/calibrated.ms', 
]

target = 'SDSSJ0901+1814'
start_freq = 35298.519 # MHz, according to 2010-2012 data
end_freq = 35298.519 + 0.5 * (256 - 1)

project_field_dict = {\
'AB1347': 'SDSS J0901+1814', 
'10C-223': 'SDSSJ0901+1814',
'12A-399': 'SDSSJ0901+1814',
'14A-371': 'J0901',
}

# split 
list_of_split_ms = []
for vis in list_of_ms:
    project_code = re.sub(r'^.*/VLA_([0-9A-Z-]+)(_[0-9A-Z-]+|)_PI_.*$', r'\1', vis)
    project_name = re.sub(r'^.*/(VLA_[0-9A-Z-_]+_PI_[^/]+)/Level_2_Calib/.*$', r'\1', vis)
    dataset_name = re.sub(r'^.*/(DataSet_[0-9]+)/calibrated/.*$', r'\1', vis)
    if project_name.find('/') >= 0:
        raise Exception('Error! Failed to parse project_name from vis string %r'%(vis))
    outputvis = 'Level_3_Split/DataSet_Merged_Ka_Band/split_%s_%s_%s.ms'%(project_name, dataset_name, target)
    print('vis', vis, 'outputvis', outputvis)
    if not os.path.isdir(os.path.dirname(outputvis)):
        os.makedirs(os.path.dirname(outputvis))
    #outputvis = vis.replace('Level_2_Calib', 'Level_3_Split').replace('calibrated/calibrated.ms', 'split_'+target+'_spw0_width4.ms')
    mstransform_params = OrderedDict()
    mstransform_params['vis'] = vis
    mstransform_params['outputvis'] = outputvis
    mstransform_params['mode'] = 'frequency'
    mstransform_params['start'] = '35298.519MHz'
    mstransform_params['width'] = '2.000MHz'
    mstransform_params['nchan'] = 64
    mstransform_params['spw'] = '*'
    mstransform_params['combinespws'] = True
    mstransform_params['regridms'] = True
    mstransform_params['nspw'] = 1
    mstransform_params['datacolumn'] = ''
    mstransform_params['field'] = target
    if not os.path.isdir(outputvis):
        # 
        tb.open(vis)
        if 'CORRECTED_DATA' in tb.colnames():
            mstransform_params['datacolumn'] = 'CORRECTED'
        else:
            mstransform_params['datacolumn'] = 'DATA'
        tb.close()
        # 
        tb.open(vis+os.sep+'SPECTRAL_WINDOW')
        spw_indicies = range(tb.nrows())
        line_spw_list = []
        for spw_index in spw_indicies:
            spw_chan_freqs = tb.getcell('CHAN_FREQ', spw_index)
            if len(spw_chan_freqs) > 1:
                min_freq = min(spw_chan_freqs)
                max_freq = max(spw_chan_freqs)
                if min_freq <= end_freq*1e6 and max_freq >= start_freq*1e6:
                    line_spw_list.append('%0d'%(spw_index))
        print('vis', vis, 'line_spw_list', line_spw_list)
        if len(line_spw_list) <= 0:
            print('Warning! The input ms data does not have a spw that contains the line!')
            continue
        mstransform_params['spw'] = ','.join(line_spw_list)
        tb.close()
        # 
        field = ''
        for key in project_field_dict:
            if vis.find(key)>=0:
                field = project_field_dict[key]
                break
        if field == '':
            raise Exception('Error! Could not find the project in project_field_dict for the vis %r. Please define this project in project_field_dict.'%(vis))
        mstransform_params['field'] = field
        # 
        print('Running mstransform(%s)'%(', '.join(['{!s}={!r}'.format(k, mstransform_params[k]) for k in mstransform_params.keys()])))
        mstransform(**mstransform_params)
        shutil.move('mstransform.last', outputvis+'.mstransform.last')
    if not os.path.isdir(outputvis):
        raise Exception('Error! Could not run split and make outputvis %r'%(outputvis))
    list_of_split_ms.append(outputvis)


# concat
vis = list_of_split_ms
concatvis = 'Level_3_Split/DataSet_Merged_Ka_Band/merged_'+target+'.ms'
if not os.path.isdir(concatvis):
    if not os.path.isdir(os.path.dirname(concatvis)):
        os.makedirs(os.path.dirname(concatvis))
    print('Running concat(vis=[%s], concatvis=%r)'%(', '.join(["'"+str(t)+"'" for t in vis]), concatvis))
    concat(vis=vis, concatvis=concatvis, freqtol='2MHz')
    shutil.move('concat.last', outputvis+'.concat.last')



