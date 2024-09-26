#!/usr/bin/bash
#SBATCH --job-name=grafi
#SBATCH --partition=day
#SBATCH --output=myjob.out
#SBATCH --error=myjob.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G

# start your logging here
echo "Starting my job"
date

# run your commands here
time /home/tadej/work/run.sh

# catch any errors here
if [ $? -ne 0 ]; then
    echo "mycommand failed"
    exit 1
fi

# end your logging here
echo "Ending my job"
date
