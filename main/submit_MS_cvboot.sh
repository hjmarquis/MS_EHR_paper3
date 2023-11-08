#!/bin/bash
#SBATCH -c 20                              # Request core
#SBATCH -N 1                               # Request node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 1-10:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=12G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=juehou@hsph.harvard.edu   # Email to which notifications will be sent
 
module load gcc/6.2.0 R/3.6.1
R CMD BATCH MS_CLIME_analysis_cvboot.R MS_CLIME_analysis_cvboot.txt