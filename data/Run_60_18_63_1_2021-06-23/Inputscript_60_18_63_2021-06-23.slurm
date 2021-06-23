#!/bin/bash -l
#SBATCH -J NUFEB_cyano1
#SBATCH -p batch
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 32
#SBATCH --cpus-per-task 1
#SBATCH --exclusive
#SBATCH --mem=100g
#SBATCH -A cnms
#SBATCH -t 10:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=None@ornl.gov

date

module purge
module load PE-gnu/3.0
#module load anaconda3

export LAMMPS=~/NUFEB-dev/src/lmp_png

ldd $LAMMPS

#run the first step to generate inputs for NUFEB
python3 ~/NUFEB-dev/nufeb-tools/GenerateAtom.py --n 10

#check if the previous run went ok, exit if not
if [ $? -ne 0 ]
then
    echo "Something went wrong in the previous step, exiting"
    exit
fi



for f in Inputscript_*.lmp
do
mpirun -np 32 $LAMMPS -in $f > ${f}.log
done

date
#do the post-processing tasks here


#cd Run_60_18_63_1
#tar -zcf VTK.tar.gz *.vtr *.vtu *.vti
#python ../../nufeb-tools/Datafedcreate.py --id None --n VTK --m ../run_60_18_63.pkl --f ./VTK.tar.gz
#python ../../nufeb-tools/Datafedcreate.py --id None --n Trajectory --m ../run_60_18_63.pkl --f ./dump.h5
#python ../../nufeb-tools/Datafedcreate.py --id None --n Log --m ../run_60_18_63.pkl --f ../NUFEB_60_18_63.log
#cd ..
