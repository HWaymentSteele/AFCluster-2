#!/bin/bash
echo $HOSTNAME
export HOME=$PWD
mkdir -p ./cache/colabfold/params/
export XDG_CACHE_HOME=./cache
export MPLCONFIGDIR=./cache

# python -m colabfold.download
cp /staging/waymentsteel/af2_weights/* ./cache/colabfold/params/
pip install scikit-learn

NAME=$1
SEQ=$2

echo ">" $NAME > ${NAME}.fasta
echo ${SEQ} >> ${NAME}.fasta

python afcluster.py --input ${NAME}.fasta --config configs/afcluster.yml