#!/bin/bash
echo $HOSTNAME
export HOME=$PWD
mkdir -p ./cache/colabfold/params/
export XDG_CACHE_HOME=./cache
export MPLCONFIGDIR=./cache

# python -m colabfold.download
cp /staging/waymentsteel/af2_weights/* ./cache/colabfold/params/
pip install scikit-learn

python afcluster.py --input INPUT.fasta --config configs/afcluster_chtc.yaml

