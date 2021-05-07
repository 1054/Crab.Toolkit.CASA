#!/bin/bash

pdoc3 --html --skip-errors -o ../../doc/html -c show_source_code=False -c show_type_annotations=True --force dzliu_combine_uvfits.py

pdoc3 --html --skip-errors -o ../../doc/html -c show_source_code=False -c show_type_annotations=True --force dzliu_clean_utils.py

