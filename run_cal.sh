#!/bin/bash
################################################################################
# UNIFIED TDDFT CALCULATION CONTROL SCRIPT
################################################################################
# This script provides centralized control for both CPU and GPU TDDFT calculations.
# Configure all settings here, then the script will automatically update and run
# the appropriate Python script (CPU or GPU).
################################################################################
# Set CUDA 12.9 library path
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:$LD_LIBRARY_PATH
# ============================================================================
# CALCULATION CONFIGURATION
# ============================================================================

# --- Molecule Selection ---
USE_XYZ=True                    # True: use XYZ file, False: use H2O test molecule
# XYZ_FILE="PTCDA_clean.xyz"      # Path to XYZ file (if USE_XYZ=True)
XYZ_FILE="H2O.xyz"      # Path to XYZ file (if USE_XYZ=True)
BASIS_SET="6-31g"               # Basis set: 6-31g, 6-31g*, def2-SVP, etc.

# --- Charge and Spin Settings ---
# CHARGE can be:
#   - Single value: CHARGE=-1
#   - Multiple values (batch): CHARGE="0 1 -1"
# For batch mode, the script will run calculations for each charge automatically
CHARGE="0 1 -1"                       # Molecular charge: 0 (neutral), +1 (cation), -1 (anion)
                                # For batch: CHARGE="0 1 -1" (space-separated)
SPIN=None                       # Spin multiplicity (2S+1): None=auto, 1=singlet, 2=doublet, 3=triplet

# --- DFT/TDDFT Settings ---
XC_FUNCTIONAL="b3lyp"           # Exchange-correlation functional: b3lyp, pbe0, cam-b3lyp, etc.
NUM_EXCITED_STATES=10           # Number of excited states to calculate
USE_TDA=False                   # True: TDA (faster), False: Full TDDFT (more accurate)

# --- Output Selection ---
# Which states to generate CUBE FILES for (0-indexed, space-separated)
# Example: "0 1 2" for first three states
STATES_TO_OUTPUT="0 1 2"

# --- Cube File Generation Options ---
GENERATE_TRANSITION_DENSITY=True    # Transition density matrix
GENERATE_EXCITED_DENSITY=True       # Excited state density
GENERATE_DENSITY_DIFFERENCE=True    # Density difference (excited - ground)
GENERATE_HOMO_LUMO=True             # HOMO and LUMO orbitals

# --- Grid Settings ---
USE_GRID_RESOLUTION=False       # True: use fixed resolution, False: use spacing
GRID_RESOLUTION_X=80            # Grid points in X (if USE_GRID_RESOLUTION=True)
GRID_RESOLUTION_Y=80            # Grid points in Y
GRID_RESOLUTION_Z=80            # Grid points in Z
BOX_MARGIN=4.0                  # Margin around molecule in Angstrom
GRID_SPACING=0.2                # Grid spacing in Angstrom

# --- NTO Analysis ---
ENABLE_NTO_ANALYSIS=True        # Generate NTO molden files
NTO_STATES="0 1 2"              # Which states for NTO analysis (0-indexed, space-separated)

# --- Transition Contribution Analysis ---
ENABLE_CONTRIBUTION_ANALYSIS=True   # Analyze orbital pair contributions
CONTRIBUTION_STATES="0 1 2"         # Which states to analyze (0-indexed, space-separated)
CONTRIBUTION_THRESHOLD=0.01         # Show contributions > 1%
TOP_N_CONTRIBUTIONS=10              # Show top N orbital pairs
GENERATE_PAIR_CUBES=True            # Generate cube files for orbital pairs
MAX_PAIRS_PER_STATE=3               # Generate cubes for top N pairs per state
PAIR_CONTRIBUTION_CUTOFF=0.05       # Only generate cubes for pairs > 5%

# --- Ground State Density and Potential ---
GENERATE_GROUND_STATE_DENSITY=True      # Generate ground state charge density cube
GENERATE_ELECTROSTATIC_POTENTIAL=True   # Generate electrostatic potential (ESP) cube

# --- Parallel Settings ---
NUM_THREADS=0                   # Number of CPU threads (0=auto-detect)

# ============================================================================
# EXECUTION CONTROL
# ============================================================================

# --- CPU or GPU Selection ---
# USE_GPU=True
USE_GPU=False

# --- Output Control ---
LOG_FILE="calculation.log"      # Log file name (auto-generated with timestamp)
AUTO_TIMESTAMP=True             # True: Add timestamp to log file name

# --- Execution Options ---
RUN_IN_BACKGROUND=False         # True: Run in background, False: Run in foreground
VERBOSE=True                    # True: Show progress, False: Quiet mode

################################################################################
# DO NOT EDIT BELOW THIS LINE (unless you know what you're doing)
################################################################################

# ============================================================================
# SCRIPT SETUP
# ============================================================================

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Determine which script to use
if [ "$USE_GPU" = True ]; then
    SCRIPT_NAME="tdm_calc_accurate_GPU.py"
    BASE_OUTPUT_DIR="output_gpu"
    CALC_TYPE="GPU"
else
    SCRIPT_NAME="tdm_calc_accurate_cpu.py"
    BASE_OUTPUT_DIR="output_cpu"
    CALC_TYPE="CPU"
fi

# Create charge-specific output directory (only for single charge mode)
# Format: output_gpu_charge-1 or output_cpu_charge0
# For batch mode, this will be set inside the loop
if [[ ! "$CHARGE" =~ [[:space:]] ]]; then
    if [ "$CHARGE" -ge 0 ]; then
        CHARGE_SUFFIX="charge${CHARGE}"
    else
        CHARGE_SUFFIX="charge${CHARGE}"  # Keeps the minus sign
    fi
    OUTPUT_DIR="${BASE_OUTPUT_DIR}_${CHARGE_SUFFIX}"
else
    # Batch mode - will be set per charge in loop
    OUTPUT_DIR="${BASE_OUTPUT_DIR}_batch"
fi

# Generate log file name with timestamp if requested (only for single charge mode)
if [[ ! "$CHARGE" =~ [[:space:]] ]]; then
    if [ "$AUTO_TIMESTAMP" = True ]; then
        TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
        LOG_FILE="${CALC_TYPE}_${CHARGE_SUFFIX}_${TIMESTAMP}.log"
    else
        LOG_FILE="${CALC_TYPE}_${CHARGE_SUFFIX}.log"
    fi
else
    # Batch mode - will be set per charge in loop
    LOG_FILE="${CALC_TYPE}_batch.log"
fi

# ============================================================================
# CONFIGURATION UPDATE FUNCTIONS
# ============================================================================

update_config() {
    local script=$1
    
    echo -e "${BLUE}Updating configuration in ${script}...${NC}"
    
    # Create a temporary Python script to update the configuration
    python3 << EOF
import re

# Read the script
with open('${script}', 'r') as f:
    content = f.read()

# Update configurations
updates = {
    'USE_XYZ': '${USE_XYZ}',
    'XYZ_FILE': "'${XYZ_FILE}'",
    'BASIS_SET': "'${BASIS_SET}'",
    'CHARGE': '${CHARGE}',
    'SPIN': '${SPIN}',
    'XC_FUNCTIONAL': "'${XC_FUNCTIONAL}'",
    'NUM_EXCITED_STATES': '${NUM_EXCITED_STATES}',
    'USE_TDA': '${USE_TDA}',
    'GENERATE_TRANSITION_DENSITY': '${GENERATE_TRANSITION_DENSITY}',
    'GENERATE_EXCITED_DENSITY': '${GENERATE_EXCITED_DENSITY}',
    'GENERATE_DENSITY_DIFFERENCE': '${GENERATE_DENSITY_DIFFERENCE}',
    'GENERATE_HOMO_LUMO': '${GENERATE_HOMO_LUMO}',
    'USE_GRID_RESOLUTION': '${USE_GRID_RESOLUTION}',
    'BOX_MARGIN': '${BOX_MARGIN}',
    'GRID_SPACING': '${GRID_SPACING}',
    'ENABLE_NTO_ANALYSIS': '${ENABLE_NTO_ANALYSIS}',
    'ENABLE_CONTRIBUTION_ANALYSIS': '${ENABLE_CONTRIBUTION_ANALYSIS}',
    'CONTRIBUTION_THRESHOLD': '${CONTRIBUTION_THRESHOLD}',
    'TOP_N_CONTRIBUTIONS': '${TOP_N_CONTRIBUTIONS}',
    'GENERATE_PAIR_CUBES': '${GENERATE_PAIR_CUBES}',
    'MAX_PAIRS_PER_STATE': '${MAX_PAIRS_PER_STATE}',
    'PAIR_CONTRIBUTION_CUTOFF': '${PAIR_CONTRIBUTION_CUTOFF}',
    'GENERATE_GROUND_STATE_DENSITY': '${GENERATE_GROUND_STATE_DENSITY}',
    'GENERATE_ELECTROSTATIC_POTENTIAL': '${GENERATE_ELECTROSTATIC_POTENTIAL}',
    'NUM_THREADS': '${NUM_THREADS}',
    'OUTPUT_DIR': "'${OUTPUT_DIR}'",
}

# Update each configuration
for key, value in updates.items():
    pattern = rf'^{key}\s*=\s*.*$'
    replacement = f'{key} = {value}'
    content = re.sub(pattern, replacement, content, flags=re.MULTILINE)

# Update list configurations (convert space-separated to Python list)
list_configs = {
    'STATES_TO_OUTPUT': '${STATES_TO_OUTPUT}',
    'NTO_STATES': '${NTO_STATES}',
    'CONTRIBUTION_STATES': '${CONTRIBUTION_STATES}',
}

for key, value in list_configs.items():
    if value.strip():
        py_list = '[' + ', '.join(value.split()) + ']'
    else:
        py_list = '[]'
    pattern = rf'^{key}\s*=\s*.*$'
    replacement = f'{key} = {py_list}'
    content = re.sub(pattern, replacement, content, flags=re.MULTILINE)

# Update GRID_RESOLUTION (special case - list of 3 values)
grid_res = f'[${GRID_RESOLUTION_X}, ${GRID_RESOLUTION_Y}, ${GRID_RESOLUTION_Z}]'
content = re.sub(r'^GRID_RESOLUTION\s*=\s*.*$', f'GRID_RESOLUTION = {grid_res}', content, flags=re.MULTILINE)

# Write back
with open('${script}', 'w') as f:
    f.write(content)

print("✓ Configuration updated successfully")
EOF
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}✗ Failed to update configuration${NC}"
        exit 1
    fi
}

# ============================================================================
# DISPLAY CONFIGURATION
# ============================================================================

display_config() {
    echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}         TDDFT CALCULATION CONFIGURATION${NC}"
    echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
    echo ""
    echo -e "${YELLOW}Execution:${NC}"
    echo "  Mode:              ${CALC_TYPE}"
    echo "  Script:            ${SCRIPT_NAME}"
    echo "  Output directory:  ${OUTPUT_DIR}"
    echo "  Log file:          ${LOG_FILE}"
    echo ""
    echo -e "${YELLOW}Molecule:${NC}"
    echo "  Use XYZ:           ${USE_XYZ}"
    if [ "$USE_XYZ" = True ]; then
        echo "  XYZ file:          ${XYZ_FILE}"
    fi
    echo "  Basis set:         ${BASIS_SET}"
    echo "  Charge:            ${CHARGE}"
    echo "  Spin:              ${SPIN}"
    echo ""
    echo -e "${YELLOW}DFT/TDDFT:${NC}"
    echo "  Functional:        ${XC_FUNCTIONAL}"
    echo "  Excited states:    ${NUM_EXCITED_STATES}"
    echo "  Use TDA:           ${USE_TDA}"
    echo ""
    echo -e "${YELLOW}Output:${NC}"
    echo "  States for cubes:  [${STATES_TO_OUTPUT}]"
    echo "  Transition dens:   ${GENERATE_TRANSITION_DENSITY}"
    echo "  Excited density:   ${GENERATE_EXCITED_DENSITY}"
    echo "  Density diff:      ${GENERATE_DENSITY_DIFFERENCE}"
    echo "  HOMO/LUMO:         ${GENERATE_HOMO_LUMO}"
    echo ""
    echo -e "${YELLOW}Analysis:${NC}"
    echo "  NTO analysis:      ${ENABLE_NTO_ANALYSIS}"
    if [ "$ENABLE_NTO_ANALYSIS" = True ]; then
        echo "    NTO states:      [${NTO_STATES}]"
    fi
    echo "  Contribution:      ${ENABLE_CONTRIBUTION_ANALYSIS}"
    if [ "$ENABLE_CONTRIBUTION_ANALYSIS" = True ]; then
        echo "    States:          [${CONTRIBUTION_STATES}]"
        echo "    Top N pairs:     ${TOP_N_CONTRIBUTIONS}"
        echo "    Pair cubes:      ${GENERATE_PAIR_CUBES}"
    fi
    echo ""
    echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
    echo ""
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

main() {
    echo ""
    display_config
    
    # Check if script exists
    if [ ! -f "${SCRIPT_NAME}" ]; then
        echo -e "${RED}✗ Error: ${SCRIPT_NAME} not found!${NC}"
        exit 1
    fi
    
    # Update configuration
    update_config "${SCRIPT_NAME}"
    echo ""
    
    # Run the calculation
    echo -e "${GREEN}Starting ${CALC_TYPE} calculation...${NC}"
    echo -e "${BLUE}Command: python3 ${SCRIPT_NAME}${NC}"
    echo -e "${BLUE}Log file: ${LOG_FILE}${NC}"
    echo ""
    
    if [ "$RUN_IN_BACKGROUND" = True ]; then
        # Run in background
        if [ "$VERBOSE" = True ]; then
            python3 "${SCRIPT_NAME}" > "${LOG_FILE}" 2>&1 &
            PID=$!
            echo -e "${GREEN}✓ Calculation started in background (PID: ${PID})${NC}"
            echo -e "${BLUE}Monitor with: tail -f ${LOG_FILE}${NC}"
        else
            python3 "${SCRIPT_NAME}" > "${LOG_FILE}" 2>&1 &
            echo -e "${GREEN}✓ Calculation started in background (PID: $!)${NC}"
        fi
    else
        # Run in foreground
        if [ "$VERBOSE" = True ]; then
            python3 "${SCRIPT_NAME}" 2>&1 | tee "${LOG_FILE}"
        else
            python3 "${SCRIPT_NAME}" > "${LOG_FILE}" 2>&1
        fi
        
        # Check exit status
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            echo ""
            echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
            echo -e "${GREEN}✓ Calculation completed successfully!${NC}"
            echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
            echo ""
            echo -e "${BLUE}Output directory: ${OUTPUT_DIR}/${NC}"
            echo -e "${BLUE}Log file: ${LOG_FILE}${NC}"
            echo ""
        else
            echo ""
            echo -e "${RED}═══════════════════════════════════════════════════════════════════${NC}"
            echo -e "${RED}✗ Calculation failed!${NC}"
            echo -e "${RED}═══════════════════════════════════════════════════════════════════${NC}"
            echo ""
            echo -e "${YELLOW}Check log file for errors: ${LOG_FILE}${NC}"
            echo ""
            exit 1
        fi
    fi
}

# ============================================================================
# BATCH MODE DETECTION AND EXECUTION
# ============================================================================

# Check if CHARGE contains multiple values (batch mode)
if [[ "$CHARGE" =~ [[:space:]] ]]; then
    # Batch mode detected
    CHARGES=($CHARGE)  # Convert to array
    TOTAL_JOBS=${#CHARGES[@]}
    
    echo ""
    echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}         BATCH MODE: ${TOTAL_JOBS} CALCULATIONS${NC}"
    echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
    echo ""
    echo -e "${YELLOW}Charges to calculate: ${CHARGE}${NC}"
    echo -e "${YELLOW}Mode: ${CALC_TYPE}${NC}"
    echo ""
    
    CURRENT_JOB=0
    FAILED_JOBS=0
    
    # Loop through each charge
    for SINGLE_CHARGE in "${CHARGES[@]}"; do
        CURRENT_JOB=$((CURRENT_JOB + 1))
        
        echo ""
        echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
        echo -e "${BLUE}JOB ${CURRENT_JOB}/${TOTAL_JOBS}: Charge = ${SINGLE_CHARGE}${NC}"
        echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
        echo ""
        
        # Temporarily set CHARGE to single value
        CHARGE=$SINGLE_CHARGE
        
        # Recalculate output directory and log file for this charge
        if [ "$SINGLE_CHARGE" -ge 0 ]; then
            CHARGE_SUFFIX="charge${SINGLE_CHARGE}"
        else
            CHARGE_SUFFIX="charge${SINGLE_CHARGE}"
        fi
        OUTPUT_DIR="${BASE_OUTPUT_DIR}_${CHARGE_SUFFIX}"
        
        if [ "$AUTO_TIMESTAMP" = True ]; then
            TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
            LOG_FILE="${CALC_TYPE}_${CHARGE_SUFFIX}_${TIMESTAMP}.log"
        else
            LOG_FILE="${CALC_TYPE}_${CHARGE_SUFFIX}.log"
        fi
        
        # Run calculation for this charge
        main
        
        # Check if successful
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✓ Job ${CURRENT_JOB}/${TOTAL_JOBS} completed successfully${NC}"
        else
            echo -e "${RED}✗ Job ${CURRENT_JOB}/${TOTAL_JOBS} failed${NC}"
            FAILED_JOBS=$((FAILED_JOBS + 1))
        fi
        
        # Small delay between jobs (if not last job)
        if [ $CURRENT_JOB -lt $TOTAL_JOBS ]; then
            sleep 1
        fi
    done
    
    # Final summary
    echo ""
    echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}         BATCH CALCULATION COMPLETE${NC}"
    echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
    echo ""
    echo -e "${YELLOW}Summary:${NC}"
    echo "  Total jobs:        ${TOTAL_JOBS}"
    echo "  Successful:        $((TOTAL_JOBS - FAILED_JOBS))"
    echo "  Failed:            ${FAILED_JOBS}"
    echo ""
    echo -e "${YELLOW}Output directories:${NC}"
    for SINGLE_CHARGE in "${CHARGES[@]}"; do
        if [ "$SINGLE_CHARGE" -ge 0 ]; then
            CHARGE_SUFFIX="charge${SINGLE_CHARGE}"
        else
            CHARGE_SUFFIX="charge${SINGLE_CHARGE}"
        fi
        echo "  Charge ${SINGLE_CHARGE}:       ${BASE_OUTPUT_DIR}_${CHARGE_SUFFIX}/"
    done
    echo ""
    
else
    # Single charge mode
    main
fi
