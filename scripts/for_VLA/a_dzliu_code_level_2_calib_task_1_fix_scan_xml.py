#!/usr/bin/env python
# 

import os, sys, glob, shutil
import xml.etree.ElementTree

input_files = glob.glob('Level_2_Calib/DataSet_*/raw/Scan.xml')

for input_file in input_files:
    
    # Check previous run
    if not os.path.isfile(input_file+'.original'):
        print('Copying "%s" to "%s"'%(input_file, input_file+'.original'))
        shutil.copy(input_file, input_file+'.original')
    
    # Open input file
    print('Reading "%s"'%(input_file))
    print('')
    tree = xml.etree.ElementTree.parse(input_file+'.original')
    
    # Loop elements
    row_nodes = tree.findall('.//row')
    has_updates = False
    for row_node in row_nodes:
        sourceName = [item.text for item in row_node.getchildren() if item.tag == 'sourceName'][0]
        scanIntent = [item.text for item in row_node.getchildren() if item.tag == 'scanIntent'][0]
        numIntent = [item.text for item in row_node.getchildren() if item.tag == 'numIntent'][0]
        calDataType = [item.text for item in row_node.getchildren() if item.tag == 'calDataType'][0]
        calibrationOnLine = [item.text for item in row_node.getchildren() if item.tag == 'calibrationOnLine'][0]
        
        if sourceName == '0137+331=3C48' and scanIntent.startswith('1 4 ') and \
            scanIntent.find('CALIBRATE_FLUX')>=0 and scanIntent.find('CALIBRATE_AMPLI')>=0 and \
            scanIntent.find('CALIBRATE_BANDPASS')>=0 and scanIntent.find('CALIBRATE_PHASE')>=0 :
            
            print('sourceName = %s'%(sourceName))
            print('scanIntent = %s'%(scanIntent))
            print('numIntent = %s'%(numIntent))
            print('calDataType = %s'%(calDataType))
            print('calibrationOnLine = %s'%(calibrationOnLine))
            
            for item in row_node.getchildren():
                if item.tag == 'numIntent':
                    item.text = '2'
                elif item.tag == 'scanIntent':
                    item.text = '1 2 CALIBRATE_BANDPASS CALIBRATE_FLUX'
                elif item.tag == 'calDataType':
                    item.text = '1 2 NONE NONE'
                elif item.tag == 'calibrationOnLine':
                    item.text = '1 2 false false'
            
            # recheck
            sourceName2 = [item.text for item in row_node.getchildren() if item.tag == 'sourceName'][0]
            scanIntent2 = [item.text for item in row_node.getchildren() if item.tag == 'scanIntent'][0]
            numIntent2 = [item.text for item in row_node.getchildren() if item.tag == 'numIntent'][0]
            calDataType2 = [item.text for item in row_node.getchildren() if item.tag == 'calDataType'][0]
            calibrationOnLine2 = [item.text for item in row_node.getchildren() if item.tag == 'calibrationOnLine'][0]
            
            if scanIntent2 == scanIntent and numIntent2 == numIntent and calDataType2 == calDataType and calibrationOnLine == calibrationOnLine2:
                print('(UNCHANGED)')
            else:
                print('(UPDATED) sourceName = %s'%(sourceName2))
                print('(UPDATED) scanIntent = %s'%(scanIntent2))
                print('(UPDATED) numIntent = %s'%(numIntent2))
                print('(UPDATED) calDataType = %s'%(calDataType2))
                print('(UPDATED) calibrationOnLine = %s'%(calibrationOnLine2))
                has_updates = True
            
            print('')
    
    
    # Write back to file
    if has_updates:
        print('Writing to "%s"'%(input_file))
        tree.write(input_file, xml_declaration = True, encoding = 'UTF-8')
    else:
        print('No update')
    
    print('')
    
    #break
