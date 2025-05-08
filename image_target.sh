#!/bin/bash --login
#SBATCH --job-name=imager
#SBATCH --partition=cosma-analyse
#SBATCH --account=do019
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=10:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

echo imager
echo $(date)

IDarray=(45589 197404 198913 202792)
#z=(0,1,2,3,4,5)

source /cosma/home/do019/dc-mcgr1/colibre_env2/bin/activate

for ID in "${IDarray[@]}"
do
    srun --exclusive -n 1 -c 1 python3 image_target.py cosma L025 m5 127 $ID &
done
wait

echo complete
echo $(date)