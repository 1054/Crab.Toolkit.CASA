#!/usr/bin/env python
# 

import os, sys, subprocess, importlib

# define function to import a package by all means
def import_package_in_casa(package_name, verbose=False):
    package = None
    try:
        aaa
        package = importlib.import_module(package_name)
    except:
        python_list = subprocess.check_output('which -a python', shell=True).split('\n')
        python_list = list(set(python_list[::-1]))[::-1]
        if verbose:
            print('python_list', python_list)
        for i in range(1, len(python_list)):
            if python_list[i] == '':
                continue
            if python_list[i].find('casa')>=0:
                continue
            command = 'PYTHONPATH="" {0} --version 2>&1'.format(python_list[i])
            if verbose:
                print('command', command)
            try:
                python_version = subprocess.check_output(command, shell=True).split('\n')[0]
            except:
                continue
            print('python_version', python_version)

            this_python_version = subprocess.check_output('python --version 2>&1', shell=True).split('\n')[0]
            print('this_python_version', this_python_version)

            if python_version.split()[1].split('.')[0] != this_python_version.split()[1].split('.')[0]:
                print('mismatched python version')
                continue

            command = 'PYTHONPATH="" {0} -c "import {1}; print({1}.__path__[0])"'.format(python_list[i], package_name)
            if verbose:
                print('command', command)
            try:
                pkg_list = subprocess.check_output(command, shell=True).split('\n')
            except:
                pkg_list = []
            if verbose:
                print('pkg_list', pkg_list)
            for pkg_path in pkg_list:
                if pkg_path == '':
                    continue
                if not os.path.exists(pkg_path):
                    continue
                pkg_dir = os.path.dirname(pkg_path)
                if verbose:
                    print('pkg_dir', pkg_dir)
                if pkg_dir not in sys.path:
                    sys.path.append(pkg_dir)
                try:
                    package = importlib.import_module(package_name)
                    if verbose:
                        print('package', package)
                except:
                    package = None
                if package is not None:
                    break
            if package is not None:
                break
    if verbose:
        print(package_name, package)
    return package

# import astropy
astropy = import_package_in_casa('astropy', verbose=False)

# 
from astropy.coordinates import SkyCoord

