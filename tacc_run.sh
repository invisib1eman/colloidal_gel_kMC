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
module load gcc/11.2.0  # If using a module system
module load gsl
module load boost
export PATH=/opt/apps/gcc/11.2.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/apps/gcc/11.2.0/lib64:$LD_LIBRARY_PATH
export PATH=/scratch/tacc/apps/gcc11_2/boost/1.86.0/include:$PATH
export LD_LIBRARY_PATH=/scratch/tacc/apps/gcc11_2/boost/1.86.0/lib:$LD_LIBRARY_PATH
export PATH=/scratch/tacc/apps/gcc11_2/gsl/2.8/include:$PATH
export LD_LIBRARY_PATH=/scratch/tacc/apps/gcc11_2/gsl/2.8/lib:$LD_LIBRARY_PATH
./colloidgel -N 8000 -L 100 -G 25 -d 0.275 -D 9mM_freeroll_step0.1_N_8000_L_100_t_300_taucrawl_100 -m 0.1 -s 300 &
./colloidgel -N 8000 -L 100 -G 25 -d 0.459 -D 3mM_freeroll_step0.1_N_8000_L_100_t_300_taucrawl_100 -m 0.1 -s 300 &
./colloidgel -N 8000 -L 100 -G 25 -d 0.6375 -D 1.7mM_freeroll_step0.1_N_8000_L_100_t_300_taucrawl_100 -m 0.1 -s 300 &
./colloidgel -N 8000 -L 100 -G 25 -d 0.919 -D 0.9mM_freeroll_step0.1_N_8000_L_100_t_300_taucrawl_100 -m 0.1 -s 300 &
./colloidgel -N 8000 -L 100 -G 25 -d 0.0925 -D 80mM_freeroll_step0.1_N_8000_L_100_t_300_taucrawl_100 -m 0.1 -s 300 &
wait
