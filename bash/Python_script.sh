#!/bin/bash -l
#$ -P mRNA velocyto
#$ -N python script
#$ -j y
#$ -m a
#$ -M yah2014@med.cornell.edu
#$ -l h_rt=4:00:00
#$ -l h_vmem=30G
#$ -q *@@red

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"  $SGE_TASK_ID"
echo "=========================================================="

file=/pbtech_mounts/homes030/yah2014/Dropbox/Public/Olivier/Projects/scRNAseq-Lymphoma/bash_python/velocyto_single.py
echo $(ls -l $file)
python $file $SGE_TASK_ID
