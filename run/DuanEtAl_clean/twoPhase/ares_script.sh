#!/bin/bash -l
## Nazwa zlecenia
#SBATCH -J FOAM
## Liczba alokowanych węzłów
#SBATCH -N 1
## Liczba zadań per węzeł (domyślnie jest to liczba alokowanych rdzeni na węźle)
#SBATCH --ntasks-per-node=24
## Ilość pamięci przypadającej na jeden rdzeń obliczeniowy (domyślnie 5GB na rdzeń)
#SBATCH --mem-per-cpu=1GB
## Maksymalny czas trwania zlecenia (format HH:MM:SS)
#SBATCH --time=71:58:00 
## Nazwa grantu do rozliczenia zużycia zasobów
#SBATCH -A plgoooopus-cpu
## Specyfikacja partycji
#SBATCH -p plgrid
## Plik ze standardowym wyjściem
#SBATCH --output="foamlog.out"
## Plik ze standardowym wyjściem błędów
#SBATCH --error="error.err"
 
 
## przejscie do katalogu z ktorego wywolany zostal sbatch
cd $SLURM_SUBMIT_DIR
srun /bin/hostname
module load openfoam/v2106-foss-2021a
source $FOAM_BASH
decomposePar > log.decomposePar -allRegions 2>&1
mpirun -np 24 multiRegionPhaseChangeFlow -parallel > log.multiRegionPhaseChangeFlow 2>&1
reconstructPar > log.reconstructPar -allRegions 2>&1
