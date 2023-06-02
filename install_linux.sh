#!/bin/bash
# This script will install SACCHARIS 2 dependencies and activate the conda environment on a linux system, and probably
# other UNIX systems as well.

# sets "unofficial bash strict mode" http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

#* check for conda and initialize environment
if command -v conda > /dev/null
then
    echo "conda found"
    conda env create -f environment.yml && . "$(conda info --base)/etc/profile.d/conda.sh" && conda activate saccharis_env \
    && conda install ./saccharis2-0.0.0dev17-py_1.tar.bz2 && echo "The SACCHARIS 2 environment is now set up and active! To run SACCHARIS first activate the environment with
     \"conda activate saccharis_env\" and then run SACCHARIS2 with \"saccharis -h\" for help"
else
    echo "Please install anaconda for your OS from https://www.anaconda.com and then run this script again"
    exit
fi
