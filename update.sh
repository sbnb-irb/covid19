#!/bin/bash
#
#

# Options for qsub
#$ -S /bin/bash
#$ -r yes
#$ -j yes
#$ -N cc_covid
#$ -wd /aloy/scratch/sbnb-www/job_covid/
#$ -m ae
#$ -M martino.bertoni@irbbarcelona.org
#$ -q ws.q
#$ -pe make 4
#$ -l mem_free=20G,h_vmem=20.2G
# End of qsub options

# Loads default environment configuration
if [[ -f $HOME/.bashrc ]]
then
  source $HOME/.bashrc
fi

OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 NUMEXPR_NUM_THREADS=4 SINGULARITYENV_PYTHONPATH=/aloy/home/mbertoni/code/chemical_checker/package SINGULARITYENV_CC_CONFIG=/aloy/home/mbertoni/cc_config.json singularity exec /aloy/data/web/sbnb_web/sbnb_web-7.59/covid19/cc_image/cc_covid19_py37.simg python /aloy/data/web/sbnb_web/sbnb_web-7.59/covid19/scripts/similarity_search.py
