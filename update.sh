#!/bin/bash
#
#

# Options for qsub
#$ -S /bin/bash
#$ -r yes
#$ -j yes
#$ -N sign4_predict_hpc.py
#$ -wd /aloy/scratch/mbertoni/tmp_jobs/tmp_bpgfsv02
#$ -pe make 4
#$ -l mem_free=10G,h_vmem=2.2G
# End of qsub options

# Loads default environment configuration
if [[ -f $HOME/.bashrc ]]
then
  source $HOME/.bashrc
fi

OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 NUMEXPR_NUM_THREADS=4 SINGULARITYENV_PYTHONPATH=/aloy/home/mbertoni/code/chemical_checker/package SINGULARITYENV_CC_CONFIG=/aloy/home/mbertoni/cc_config.json singularity exec /aloy/data/web/sbnb_web/sbnb_web-7.59/covid19/cc_image/cc_covid19_py37.simg python /aloy/data/web/sbnb_web/sbnb_web-7.59/covid19/scripts/similarity_search.py