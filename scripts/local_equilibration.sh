#!/bin/bash
# Local equilibration script for OX1R simulation on RTX 3060 GPU.
# Execute this from the project root.

set -e # Exit on error

# 1. GROMACS Environment Activation
# Source your high-performance CUDA build
if [ -f "$HOME/gromacs-cuda/bin/GMXRC" ]; then
    echo "Activating High-Performance CUDA GROMACS..."
    source "$HOME/gromacs-cuda/bin/GMXRC"
elif command -v gmx &> /dev/null; then
    echo "Using system/environment GROMACS: $(which gmx)"
else
    echo "Error: GROMACS not found. Please ensure it is installed and in your PATH."
    exit 1
fi

# --- SETTINGS ---
USE_GPU=true
CPU_CORES=$(nproc)
# ----------------

echo "Mode: GPU Accelerated (CUDA)"
echo "Starting local equilibration on $(hostname)..."

# Canonical project locations
INPUT_COORDS_DIR="system/coordinates"
TOPOLOGY_DIR="system/topology"
MDP_DIR="system/mdp"
RUN_DIR="runs/equilibration"

# 2. Run Stepwise Equilibration (0 to 6)
prev="input"

for step in 0_minimization 1_equilibration 2_equilibration 3_equilibration 4_equilibration 5_equilibration 6_equilibration
do
    echo "================================================================================"
    echo " PHASE: ${step}"
    echo "================================================================================"
    
    if [ "$prev" = "input" ]; then
        prev_gro="${INPUT_COORDS_DIR}/${prev}.gro"
    else
        prev_gro="${RUN_DIR}/${prev}.gro"
    fi

    if [ ! -f "${prev_gro}" ]; then
        echo "Error: Input file ${prev_gro} not found."
        exit 1
    fi

    gmx grompp -f "${MDP_DIR}/${step}.mdp" \
               -c "${prev_gro}" \
               -p "${TOPOLOGY_DIR}/topol.top" \
               -n "${TOPOLOGY_DIR}/index.ndx" \
               -o "${RUN_DIR}/${step}.tpr" \
               -r "${prev_gro}" \
               -maxwarn 2
    
    # Run the MD engine
    if [[ "$step" == *"minimization"* ]]; then
        echo "Running Minimization (GPU for non-bonded)..."
        gmx mdrun -v -deffnm "${RUN_DIR}/${step}" -nb gpu -nt $CPU_CORES
    else
        echo "Running Equilibration (GPU for Non-bonded and PME)..."
        # Removed -update gpu because the system contains virtual sites
        gmx mdrun -v -deffnm "${RUN_DIR}/${step}" \
                  -nb gpu \
                  -pme gpu \
                  -nt $CPU_CORES
    fi
    
    prev=${step}
    echo "${step} completed successfully."
done

echo "================================================================================"
echo "All equilibration steps finished!"
echo "Final structure: ${RUN_DIR}/6_equilibration.gro"
echo "Ready to transfer to cluster."
echo "================================================================================"
