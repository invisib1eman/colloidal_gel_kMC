#!/bin/sh
#----------------------------------------------------

#SBATCH -J my0.01_10.0_2.0_.49
#SBATCH -o my.o%j       # Name of stdout output file
#SBATCH -e my.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 128              # Total # of mpi tasks
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A OTH21040       # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=dwqian@utexas.edu

# Any other commands must follow all #SBATCH directives...
cd $SLURM_SUBMIT_DIR
module load gsl
module load boost
./colloidgel -N 8000 -L 100 -G 25 -d 0.275 -D 9mM_freeroll_step0.1_N_8000_L_100 -m 0.1 &
./colloidgel -N 8000 -L 100 -G 25 -d 0.459 -D 3mM_freeroll_step0.1_N_8000_L_100 -m 0.1 &
./colloidgel -N 8000 -L 100 -G 25 -d 0.6375 -D 1.7mM_freeroll_step0.1_N_8000_L_100 -m 0.1 &
./colloidgel -N 8000 -L 100 -G 25 -d 0.919 -D 0.9mM_freeroll_step0.1_N_8000_L_100 -m 0.1 &
./colloidgel -N 8000 -L 100 -G 25 -d 0.0925 -D 80mM_freeroll_step0.1_N_8000_L_100 -m 0.1 &
wait
