#!/bin/bash
# Local production script for the prepared MD system.
# Execute this from the project root.

set -e

if [ -f "/home/jallow/gromacs-2026.1/build/scripts/GMXRC" ]; then
    source "/home/jallow/gromacs-2026.1/build/scripts/GMXRC"
elif command -v gmx &> /dev/null; then
    true
else
    echo "Error: GROMACS not found. Please ensure it is installed and in your PATH."
    exit 1
fi

CPU_CORES=$(nproc)
RUN_DIR="runs/production"
TPR_FILE="${RUN_DIR}/production.tpr"

if [ ! -f "${TPR_FILE}" ]; then
    echo "Error: ${TPR_FILE} not found."
    echo "Prepare it with:"
    echo "gmx grompp -f system/mdp/production.mdp -c runs/equilibration/6_equilibration.gro -t runs/equilibration/6_equilibration.cpt -p system/topology/topol.top -n system/topology/index.ndx -o runs/production/production.tpr"
    exit 1
fi

gmx mdrun -v -deffnm "${RUN_DIR}/production" \
          -nb gpu \
          -pme gpu \
          -nt "${CPU_CORES}"
