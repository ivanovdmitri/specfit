#!/bin/bash

unset SPECFIT
if [ "x${BASH_ARGV[0]}" = "x" ]; then
    if [ ! -f ./specfit_env.sh ]; then
        echo ERROR: must "cd where/specfit/is" before calling "./specfit_env.sh" for this version of bash!
        SPECFIT=; export SPECFIT
        return
    fi
    SPECFIT="$PWD"; export SPECFIT
else
    # get param to "."
    THIS=$(dirname ${BASH_ARGV[0]})
    SPECFIT=$(cd ${THIS};pwd); export SPECFIT
fi

LD_LIBRARY_PATH=$SPECFIT/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

PYTHONPATH=$SPECFIT:$PYTHONPATH
export PYTHONPATH

