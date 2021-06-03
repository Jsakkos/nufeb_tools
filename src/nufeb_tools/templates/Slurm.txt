#!/bin/bash -l
#SBATCH -J ${job}
#SBATCH -p batch
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 32
#SBATCH --cpus-per-task 1
#SBATCH --exclusive
#SBATCH --mem=100g
#SBATCH -A cnms
#SBATCH -t 10:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${USER}@ornl.gov

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


#cd Run_${n_cyanos}_${n_ecw}_${SucPct}_${Replicates}
#tar -zcf VTK.tar.gz *.vtr *.vtu *.vti
#python ../../nufeb-tools/Datafedcreate.py --id ${id} --n VTK --m ../run_${n_cyanos}_${n_ecw}_${SucPct}.pkl --f ./VTK.tar.gz
#python ../../nufeb-tools/Datafedcreate.py --id ${id} --n Trajectory --m ../run_${n_cyanos}_${n_ecw}_${SucPct}.pkl --f ./dump.h5
#python ../../nufeb-tools/Datafedcreate.py --id ${id} --n Log --m ../run_${n_cyanos}_${n_ecw}_${SucPct}.pkl --f ../NUFEB_${n_cyanos}_${n_ecw}_${SucPct}.log
#cd ..
