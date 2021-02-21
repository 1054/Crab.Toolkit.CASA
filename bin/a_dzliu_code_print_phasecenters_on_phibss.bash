#!/bin/bash

source ~/Software/CASA/SETUP.bash 5.4.0

casa --nogui --log2term -c "execfile('a_dzliu_code_print_phasecenters_in_casa.py')"


