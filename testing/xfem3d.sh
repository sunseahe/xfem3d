#!/bin/bash
#PBS -m e
#PBS -j oe
#PBS -r n
#PBS -q cn1cpu16
#PBS -l nodes=1:ppn=16
#PBS -l select=1
cd ${PBS_O_WORKDIR}
../xfem -i unit_brick_r3 -c 16
