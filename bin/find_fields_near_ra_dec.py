#!/usr/bin/env python
# 
# must have casa in PATH
# 

import os, sys, re, math, subprocess

casa_executables = []
if len(casa_executables) == 0:
    casa_executables = subprocess.check_output('which casa 2>/dev/null; echo ""; exit 0', shell=True).split('\n')
    casa_executables = [t for t in casa_executables if t != '']
if len(casa_executables) == 0:
    casa_executables = subprocess.check_output('alias | grep casa | tac', shell=True).split('\n')
    print('casa_executables', casa_executables)
    casa_executables = [re.sub(r'^ alias casa*=\'(.*)\'', r'\1', t) for t in casa_executables if t != '']
    print('casa_executables', casa_executables)
if len(casa_executables) == 0:
    raise Exception('Error! Could not find casa executable in PATH!')
casa_executable = casa_executables[0]


# define function to check casa enviornment
def check_casa():
    if sys.argv[0].endswith('start_casa.py') or sys.argv[0].endswith('casashell/__main__.py'):
        return True
    else:
        return False


# define function to convert RA Dec from sexagesimal to decimal
def convert_sexagesimal_to_decimal(ra_str, dec_str, verbose=False):
    ra_regexmatch = re.match('^([0-9+-]+)h([0-9]+)m([0-9.]+)s', ra_str)
    ra_hour = None
    if ra_regexmatch:
        ra_hour, ra_min, ra_sec = ra_regexmatch.groups()
    else:
        ra_regexmatch = re.match('^([0-9+-]+):([0-9]+):([0-9.]+)', ra_str)
        if ra_regexmatch:
            ra_hour, ra_min, ra_sec = ra_regexmatch.groups()
    
    ra_degrees = None
    if ra_hour is not None:
        if ra_hour[0] == '-':
            ra_degrees = -( math.fabs(float(ra_hour)) + float(ra_min)/60.0 + float(ra_sec)/3600.0 ) * 15.0
        else:
            ra_degrees = ( math.fabs(float(ra_hour)) + float(ra_min)/60.0 + float(ra_sec)/3600.0 ) * 15.0

    dec_deg = None
    dec_regexmatch = re.match('^([0-9+-]+)d([0-9]+)m([0-9.]+)s', dec_str)
    if dec_regexmatch:
        dec_deg, dec_min, dec_sec = dec_regexmatch.groups()
    else:
        dec_regexmatch = re.match('^([0-9+-]+):([0-9]+):([0-9.]+)', dec_str)
        if dec_regexmatch:
            dec_deg, dec_min, dec_sec = dec_regexmatch.groups()
        else:
            dec_regexmatch = re.match(r'^([0-9+-]+)\.([0-9]+)\.([0-9.]+)', dec_str)
            if dec_regexmatch:
                dec_deg, dec_min, dec_sec = dec_regexmatch.groups()
    
    dec_degrees = None
    if dec_deg is not None:
        if dec_deg[0] == '-':
            dec_degrees = -( math.fabs(float(dec_deg)) + float(dec_min)/60.0 + float(dec_sec)/3600.0 )
        else:
            dec_degrees = ( math.fabs(float(dec_deg)) + float(dec_min)/60.0 + float(dec_sec)/3600.0 )

    if verbose:
        print('{} {} -> {} {}'.format(ra_str, dec_str, ra_degrees, dec_degrees))

    return ra_degrees, dec_degrees

# 
#ra, dec = convert_sexagesimal_to_decimal('17:45:40.307861', '-29.00.34.42600', verbose=True)


# 
if check_casa():
    #subprocess.check_output("PYTHONPATH=\"\" %s "%(sys.argv[0]), shell=True)
    #sys.exit()
    raise Exception("Please run this script outside CASA!")


# 
#if os.path.exists('/software/casa/analysis_scripts'):
#    sys.path.append('/software/casa/analysis_scripts')
#import analysisUtils as aU
#vm = aU.ValueMapping(vis)
#vm.uniqueFields
#mytb = aU.createCasaTool(aU.tbtool)
#mytb.open(vis+'/FIELD')
#mytb.close()


# Read User Input
if len(sys.argv) != 3+1:
    print('Usage: find_fields_near_ra_dec.py <vis> <ra> <dec>')
    sys.exit()

vis = sys.argv[1]
ura = sys.argv[2]
udec = sys.argv[3]


# Write CASA script
if not os.path.exists(vis+'.fields.ra.dec.txt'):
    script_content = "import os, math\n"
    script_content += "vis = '%s'\n"%(vis)
    script_content += """\
    tb.open(vis+'/FIELD')
    name_list = tb.getcol('NAME')
    ra_list, dec_list = tb.getcol('PHASE_DIR')
    if ra_list.shape[0] == 1:
        ra_list = ra_list[0]
    if dec_list.shape[0] == 1:
        dec_list = dec_list[0]
    ra_list = ra_list/math.pi*180.0+(ra_list<0)*360.0
    dec_list = dec_list/math.pi*180.0
    tb.close()
    with open(vis+'.fields.ra.dec.txt', 'w') as fp:
        fp.write("# %-16s %15s %15s\\n"%("Name", "RA", "DEC"))
        for i in range(len(name_list)):
            fp.write("%-18s %15.7f %15.7f\\n"%(name_list[i], ra_list[i], dec_list[i]))
    print("Output to %s"%(vis+'.fields.ra.dec.txt'))
    """

    with open(vis+'.fields.ra.dec.script', 'w') as fp:
        fp.write(script_content+"\n")

    subprocess.check_output("%s --nogui --nologger -c 'execfile(\"%s\")'"%(casa_executable, vis+'.fields.ra.dec.script'), 
                            shell=True)

name_list = []
ra_list = []
dec_list = []
with open(vis+'.fields.ra.dec.txt', 'r') as fp:
    for line in fp.readlines():
        if line.startswith('#'):
            continue
        linesplit = line.strip().split()
        if len(linesplit) == 3:
            name_list.append(linesplit[0])
            ra_list.append(float(linesplit[1]))
            dec_list.append(float(linesplit[2]))


# get ra dec in degrees
try:
    float(ura)
    float(udec)
    ra, dec = float(ura), float(udec)
except:
    ra, dec = convert_sexagesimal_to_decimal(ura, udec, verbose=True)


# 
idx_list = []
sep_list = []
for i in range(len(name_list)):
    #print('ra dec list', ra_list[i], dec_list[i], 'ra dec', ra, dec)
    sep_arcsec = math.sqrt((ra_list[i]-ra)**2 + (dec_list[i]-dec)**2) * 3600.0
    sep_list.append(sep_arcsec)
    idx_list.append(i)

_, sorted_idx_list = zip(*sorted(zip(sep_list, idx_list)))

print("# %-16s %15s %15s %15s"%("Name", "RA", "DEC", "SEP[arcsec]"))
for i in sorted_idx_list:
    print("%-18s %15.7f %15.7f %15.3f"%(name_list[i], ra_list[i], dec_list[i], sep_list[i]))


