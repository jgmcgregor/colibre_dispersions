#!/bin/bash --login
#SBATCH --job-name=multiselect
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

snapArray=(0127 0092 0076 0064 0056 0048)
#z=(0,1,2,3,4,5)

source /cosma/apps/do019/dc-mcgr1/colibre_env/bin/activate

for snap in "${snapArray[@]}"
do
    srun --nodes=1 --ntasks=1 python3 -Wignore selection.py cosma L025 m5 $snap zero &
done
wait

echo complete
echo $(date)
