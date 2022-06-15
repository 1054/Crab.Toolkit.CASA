# 
# -*- coding: utf-8 -*-
# https://stackoverflow.com/questions/436198/what-is-an-alternative-to-execfile-in-python-3
# 

def execfile(filepath, globals=None, locals=None):
    if globals is None:
        globals = {}
    globals.update({
        "__file__": filepath,
        "__name__": "__main__",
    })
    with open(filepath, 'rb') as file:
        exec(compile(file.read(), filepath, 'exec'), globals, locals)

# execute the file
#execfile("/path/to/somefile.py")

