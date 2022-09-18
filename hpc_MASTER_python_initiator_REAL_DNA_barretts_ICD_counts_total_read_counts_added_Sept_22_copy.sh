
#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=16G
#$ -l h_rt=24:0:0
#$ -cwd
#$ -j y

export OMP_NUM_THREADS=$NSLOTS

module load python
source /data/home/hfx381/presto/bin/activate

## EDIT - Ensure the python script has EXECUTE permissions byt typing:
## $chmod u+x file

#PATH to the python functions

pythonfiles="/data/BCI-ESCS/richard/processing_files/python_scripts"

## EDIT the name of the python script for the databatch txt files that are to be processed

$pythonfiles/hpc_master_combined_script_for_nan_correction_REAL_DNA_barretts_ICD_counts_total_read_counts_added_Sept_22.py 

deactivate

module unload python

exit
