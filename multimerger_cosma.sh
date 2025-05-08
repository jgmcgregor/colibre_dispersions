#!/bin/bash --login
#SBATCH --job-name=multimerger
#SBATCH --partition=cosma-analyse
#SBATCH --account=do019
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1g
#SBATCH --time=10:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

echo multiselect
echo $(date)

snapArray=(127 92 76 64 56 48)
#z=(0,1,2,3,4,5)

source /cosma/home/do019/dc-mcgr1/colibre_env2/bin/activate

for snap in "${snapArray[@]}"
do
    srun --exclusive -n 1 -c 1 python3 mergercheck.py cosma L025 m5 $snap &
done
wait

echo complete
echo $(date)
