#!/bin/bash
# 

# source certain CASA version if on certain machine
if [[ -f "$HOME/Software/CASA/SETUP.bash" ]]; then
    source "$HOME/Software/CASA/SETUP.bash" 5.7.2
fi

# check CASA
if [[ $(type casa 2>/dev/null | wc -l) -eq 0 ]]; then
    echo "Error! Could not find casa in \$PATH!"
    exit 255
fi

# read user input and export global variable
usage() {
    echo "Usage:"
    echo "    ./a_dzliu_code_run_casa_tclean_wide_line_cube.bash -vis your_input_vis.ms -out your_output_prefix"
    echo "Notes:"
    echo "    Optional argument -target or -field can specify a field name in the input ms data."
}
CASA_COMMANDLINE_VIS=""
CASA_COMMANDLINE_OUT=""
CASA_COMMANDLINE_FIELD=""
iarg=1
argstr=""
argmode="none"
while [[ $iarg -le $# ]]; do
    argstr=$(echo "${!iarg}" | perl -p -e 's/^[-]+/-/g' | tr [:upper:] [:lower:])
    case "$argstr" in
        -vis)
            argmode="CASA_COMMANDLINE_VIS"; iarg=$((iarg+1))
            ;;
        -out)
            argmode="CASA_COMMANDLINE_OUT"; iarg=$((iarg+1))
            ;;
        -target|-field)
            argmode="CASA_COMMANDLINE_FIELD"; iarg=$((iarg+1))
            ;;
        -*)
            if [[ $(echo "$argstr" | perl -p -e 's/^-[0-9.]([0-9eE.]*)$/ISNUMBER/g') != "ISNUMBER" ]]; then
                argmode="none"
                echo "Warning! Input argument ${!iarg} is not recognized!"
            fi
            ;;
    esac
    if [[ $iarg -gt $# ]]; then
        break
    fi
    if [[ "$argmode" != "none" ]]; then
        ## for array
        #if [[ "${!argmode}"x != ""x ]]; then
        #    varcontent=$(eval printf \"\\\"%q\\\" \" \${${argmode}[@]})
        #else
        #    varcontent=""
        #fi
        #echo "$argmode=(${varcontent}\"${!iarg}\")"
        #eval "$argmode=(${varcontent}\"${!iarg}\")"
        echo "$argmode=\"${!iarg}\""
        eval "$argmode=\"${!iarg}\""
    fi
    iarg=$((iarg+1))
done
#echo "CASA_COMMANDLINE_VIS = ${CASA_COMMANDLINE_VIS[@]} (${#CASA_COMMANDLINE_VIS[@]})"
#echo "CASA_COMMANDLINE_FIELD = ${CASA_COMMANDLINE_FIELD[@]} (${#CASA_COMMANDLINE_FIELD[@]})"
#echo "CASA_COMMANDLINE_OUT = ${CASA_COMMANDLINE_OUT[@]} (${#CASA_COMMANDLINE_OUT[@]})"
#if [[ ${#CASA_COMMANDLINE_VIS[@]} -eq 0 ]] || [[ ${#CASA_COMMANDLINE_OUT[@]} -eq 0 ]]; then
#    usage
#    exit 255
#fi
if [[ "$CASA_COMMANDLINE_VIS"x == ""x ]] || [[ "$CASA_COMMANDLINE_OUT"x == x ]]; then
    usage
    exit 255
fi
export CASA_COMMANDLINE_VIS
export CASA_COMMANDLINE_FIELD
export CASA_COMMANDLINE_OUT
#exit

# run python script in CASA
script_dir=$(dirname $(perl -MCwd -e 'print Cwd::abs_path shift' "${BASH_SOURCE[0]}"))
script_name=$(basename "${BASH_SOURCE[0]}" | perl -p -e 's/\.bash$//g')
echo "casa --nogui --log2term -c \"exec(compile(open(\\\"${script_dir}/${script_name}.py\\\").read(),\\\"${script_dir}/${script_name}.py\\\",\\\"exec\\\"),globals(),locals())\""
casa --nogui --log2term -c "exec(compile(open(\"${script_dir}/${script_name}.py\").read(),\"${script_dir}/${script_name}.py\",\"exec\"),globals(),locals())"

#locals().update(dict(input_vis='aaa',output_prefix='bbb'))


