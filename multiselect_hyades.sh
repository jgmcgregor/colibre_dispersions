#!/bin/bash --login
#SBATCH --job-name=multiselect
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8g
#SBATCH --time=10:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

echo multiselect
echo $(date)

snapArray=(0123 0088 0072 0064 0056 0048)
#z=(0,1,2,3,4,5)

source /home/jmcgregor/envs/colibre_env/bin/activate

for snap in "${snapArray[@]}"
do
    srun --nodes=1 --ntasks=1 python3 -Wignore selection.py hyades L0025 N0752 $snap zero &
done
wait

echo complete
echo $(date)
