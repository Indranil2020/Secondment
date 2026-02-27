#!/bin/bash
################################################################################
# UNIFIED DFT/TDDFT CALCULATION CONTROL SCRIPT
################################################################################
# This script provides centralized control for both CPU and GPU calculations.
# Configure all settings here, then the script will automatically update and run
# the appropriate Python script (CPU or GPU).
#
# WORKFLOW: Input Structure → [Optimization] → [DFT] → [TDDFT] → Output
#           Each stage can be enabled/disabled independently
################################################################################

# Set CUDA 12.9 library path (for GPU calculations)
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:$LD_LIBRARY_PATH

# ============================================================================
#                    SECTION 1: CALCULATION STAGE CONTROL
# ============================================================================
# Control which calculation stages to run. Valid combinations:
#   Opt only:      OPTIMISE_GEOMETRY=True,  ENABLE_DFT=False, ENABLE_TDDFT=False
#   DFT only:      OPTIMISE_GEOMETRY=False, ENABLE_DFT=True,  ENABLE_TDDFT=False
#   Opt + DFT:     OPTIMISE_GEOMETRY=True,  ENABLE_DFT=True,  ENABLE_TDDFT=False
#   DFT + TDDFT:   OPTIMISE_GEOMETRY=False, ENABLE_DFT=True,  ENABLE_TDDFT=True
#   Full workflow: OPTIMISE_GEOMETRY=True,  ENABLE_DFT=True,  ENABLE_TDDFT=True
# Note: TDDFT requires DFT (auto-enabled if ENABLE_TDDFT=True)
# ----------------------------------------------------------------------------
ENABLE_DFT=True                 # Run ground state DFT calculation
ENABLE_TDDFT=True               # Run TDDFT excited state calculation

# ============================================================================
#                    SECTION 2: INPUT STRUCTURE
# ============================================================================
USE_XYZ=True                    # True: use XYZ file, False: use built-in H2O test
# XYZ_FILE="PTCDA_clean.xyz"      # Path to XYZ file (relative or absolute)
XYZ_FILE="opt/charge0/wb97x_d3bj/optimised_structure.xyz"     # PTCDA Neutral optimised structure
# XYZ_FILE="emission/optimised_structure_wb97x_d3bj_6_31g__gpu_charge0/optimised_structure.xyz"
# XYZ_FILE="emission/charge-1_s1_opt/optimised_structure.xyz"     # PTCDA Anion optimised structure
# XYZ_FILE="/home/indranil/Documents/Secondment/test/PTCDA/anion_cal/opt/charge-1/wb97x_d3bj/ptcda_anion_dimer_T_from_molecules_ini.xyz"
# XYZ_FILE="azulene.xyz"
# XYZ_FILE="azulene_-1_opt.xyz"

# ============================================================================
#                    SECTION 3: ELECTRONIC STRUCTURE
# ============================================================================
# --- Charge and Spin ---
# CHARGE="-2"                    # Molecular charge: 0, +1, -1, or batch "0 1 -1"
# SPIN=3                          # Spin multiplicity: 3=triplet for dianion (2 unpaired e-)
CHARGE="0"
SPIN=None

# --- Constrained DFT (cDFT) ---
# Use this for multi-fragment systems to enforce a target charge on a fragment (Monomer A).
# Notes on functional/basis choices for charge-transfer and cDFT-style constraints:
# - Best practice (often): range-separated hybrids like wb97x, cam-b3lyp for CT states.
# - Basis: include diffuse functions for anions/CT: def2-SVPD/def2-TZVPPD/def2-QZVPPD.
ENABLE_CDFT=False
# 0-indexed atom indices belonging to Monomer A (space-separated). Example: "0 1 2 3 4"
MONOMER_A_ATOMS=""
# Target charge on Monomer A. Example: -1 enforces one extra electron on fragment A.
TARGET_CHARGE_A=0.0
# Secant initial guesses for constraint multiplier Vc
CDFT_VC_X0=0.0
CDFT_VC_X1=0.1
# Convergence and iteration limits
CDFT_CHARGE_TOL=1e-4
CDFT_MAX_ITER=25

# --- Basis Set ---
BASIS_SET="def2-TZVPPD"              # Options: 6-31g, 6-31g*, 6-311g**, def2-SVP, def2-TZVP, def2-SVPD, def2-QZVPPD 

# --- DFT Functional ---
# Supported dispersion: -d3bj, -d3zero, -d3bjm, -d3zerom, -d3op, -d4
# NOT supported: wb97x-d, wb97x-d3 (use wb97x-d3bj instead)
XC_FUNCTIONAL="wb97x-d3bj"      # Examples: b3lyp, pbe0, cam-b3lyp, wb97x-d3bj

# ============================================================================
#                    SECTION 4: GEOMETRY OPTIMIZATION
# ============================================================================
OPTIMISE_GEOMETRY=True           # True: optimize geometry, False: use input as-is
OPT_CYCLES=5                    # Number of optimization cycles (output N → input N+1)
OPT_MAX_STEPS=150               # Max steps per optimization cycle
OPT_CONV_PARAMS="tight"         # Convergence: "tight", "normal", "loose"
                                #   tight:  energy 1e-6, gradient 3e-4
                                #   normal: energy 1e-5, gradient 1e-3
                                #   loose:  energy 1e-4, gradient 3e-3

# ============================================================================
#                    SECTION 5: TDDFT SETTINGS
# ============================================================================
NUM_EXCITED_STATES=6          # Number of excited states to calculate (0 = skip TDDFT)
USE_TDA=False                   # True: TDA (faster), False: Full TDDFT (more accurate)

# --- Emission Calculation ---
# Emission = fluorescence from excited state minimum geometry
# Requires excited state geometry optimization + TDDFT at that geometry
ENABLE_EMISSION=False           # Calculate emission energy (requires ENABLE_TDDFT=True)
EMISSION_STATE=2                # Which state to optimize (0=S1, 1=S2, etc.)
EMISSION_OPT_MAX_STEPS=200      # Max steps for excited state optimization
EMISSION_OPT_CONV="normal"       # Convergence: "tight", "normal", "loose"

# ============================================================================
#                    SECTION 6: OUTPUT GENERATION
# ============================================================================
# --- Ground State Outputs ---
GENERATE_GROUND_STATE_DENSITY=True      # Ground state electron density cube
GENERATE_ELECTROSTATIC_POTENTIAL=True   # Electrostatic potential (ESP) cube
GENERATE_DEFORMATION_DENSITY=True       # Deformation density (SCF - Promolecule)
GENERATE_HOMO_LUMO=True                 # HOMO/LUMO orbital cubes

# --- Excited State Outputs ---
STATES_TO_OUTPUT="0 1 2 3 4 5"                # States for cube files (0-indexed, space-separated)
GENERATE_TRANSITION_DENSITY=True        # Transition density matrix cubes
GENERATE_EXCITED_DENSITY=True           # Excited state density cubes
GENERATE_DENSITY_DIFFERENCE=True        # Density difference cubes (excited - ground)

# --- NTO Analysis ---
ENABLE_NTO_ANALYSIS=True                # Generate NTO molden files
NTO_STATES="0 1 2 3 4 5"                      # States for NTO analysis (0-indexed)

# --- Transition Contribution Analysis ---
ENABLE_CONTRIBUTION_ANALYSIS=True       # Analyze orbital pair contributions
CONTRIBUTION_STATES="0 1 2 3 4 5"             # States to analyze
CONTRIBUTION_THRESHOLD=0.01             # Show contributions > 1%
TOP_N_CONTRIBUTIONS=5                  # Show top N orbital pairs it was 10
GENERATE_PAIR_CUBES=True                # Generate cube files for orbital pairs
MAX_PAIRS_PER_STATE=3                   # Max cubes per state
PAIR_CONTRIBUTION_CUTOFF=0.05           # Only pairs > 5%

# ============================================================================
#                    SECTION 7: GRID SETTINGS (for cube files)
# ============================================================================
USE_GRID_RESOLUTION=False       # True: fixed resolution, False: use spacing
GRID_RESOLUTION_X=80            # Grid points in X (if USE_GRID_RESOLUTION=True)
GRID_RESOLUTION_Y=80            # Grid points in Y
GRID_RESOLUTION_Z=80            # Grid points in Z
BOX_MARGIN=15                  # Margin around molecule (Angstrom)
GRID_SPACING=0.15                # Grid spacing (Angstrom) - fine grid for accurate integration

# ============================================================================
#                    SECTION 8: EXECUTION SETTINGS
# ============================================================================
# --- Hardware Selection ---
USE_GPU=True                    # True: GPU (faster), False: CPU

# --- Parallel Settings ---
NUM_THREADS=0                   # CPU threads (0=auto-detect)

# --- Verbosity ---
VERBOSE_LEVEL=4                 # 0: minimal, 1: normal, 2: detailed, 3: debug, 4: max

# --- Output Control ---
LOG_FILE="calculation.log"      # Log file name
AUTO_TIMESTAMP=True             # Add timestamp to log file name
RUN_IN_BACKGROUND=False         # Run in background
VERBOSE=True                    # Show progress

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

# Extract molecule name from XYZ_FILE (without .xyz extension)
MOLECULE_NAME=$(basename "$XYZ_FILE" .xyz)

# Clean up XC functional name for directory (replace special chars)
XC_CLEAN=$(echo "$XC_FUNCTIONAL" | tr '[:upper:]' '[:lower:]' | sed 's/[^a-z0-9]/_/g')
# Clean up basis set name for directory
BASIS_CLEAN=$(echo "$BASIS_SET" | tr '[:upper:]' '[:lower:]' | sed 's/[^a-z0-9]/_/g')

# Determine which script to use
if [ "$USE_GPU" = True ]; then
    SCRIPT_NAME="tdm_calc_accurate_GPU.py"
    BASE_OUTPUT_DIR="${MOLECULE_NAME}_${XC_CLEAN}_${BASIS_CLEAN}_gpu"
    CALC_TYPE="GPU"
else
    SCRIPT_NAME="tdm_calc_accurate_cpu.py"
    BASE_OUTPUT_DIR="${MOLECULE_NAME}_${XC_CLEAN}_${BASIS_CLEAN}_cpu"
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
    # Section 1: Calculation Stage Control
    'ENABLE_DFT': '${ENABLE_DFT}',
    'ENABLE_TDDFT': '${ENABLE_TDDFT}',
    # Section 2: Input Structure
    'USE_XYZ': '${USE_XYZ}',
    'XYZ_FILE': "'${XYZ_FILE}'",
    # Section 3: Electronic Structure
    'BASIS_SET': "'${BASIS_SET}'",
    'CHARGE': '${CHARGE}',
    'SPIN': '${SPIN}',
    'XC_FUNCTIONAL': "'${XC_FUNCTIONAL}'",
    # Constrained DFT
    'ENABLE_CDFT': '${ENABLE_CDFT}',
    'TARGET_CHARGE_A': '${TARGET_CHARGE_A}',
    'CDFT_VC_X0': '${CDFT_VC_X0}',
    'CDFT_VC_X1': '${CDFT_VC_X1}',
    'CDFT_CHARGE_TOL': '${CDFT_CHARGE_TOL}',
    'CDFT_MAX_ITER': '${CDFT_MAX_ITER}',
    # Section 4: Geometry Optimization
    'OPTIMISE_GEOMETRY': '${OPTIMISE_GEOMETRY}',
    'OPT_CYCLES': '${OPT_CYCLES}',
    'OPT_MAX_STEPS': '${OPT_MAX_STEPS}',
    'OPT_CONV_PARAMS': "'${OPT_CONV_PARAMS}'",
    # Section 5: TDDFT Settings
    'NUM_EXCITED_STATES': '${NUM_EXCITED_STATES}',
    'USE_TDA': '${USE_TDA}',
    # Section 5b: Emission Calculation
    'ENABLE_EMISSION': '${ENABLE_EMISSION}',
    'EMISSION_STATE': '${EMISSION_STATE}',
    'EMISSION_OPT_MAX_STEPS': '${EMISSION_OPT_MAX_STEPS}',
    'EMISSION_OPT_CONV': "'${EMISSION_OPT_CONV}'",
    # Section 6: Output Generation
    'GENERATE_TRANSITION_DENSITY': '${GENERATE_TRANSITION_DENSITY}',
    'GENERATE_EXCITED_DENSITY': '${GENERATE_EXCITED_DENSITY}',
    'GENERATE_DENSITY_DIFFERENCE': '${GENERATE_DENSITY_DIFFERENCE}',
    'GENERATE_HOMO_LUMO': '${GENERATE_HOMO_LUMO}',
    'ENABLE_NTO_ANALYSIS': '${ENABLE_NTO_ANALYSIS}',
    'ENABLE_CONTRIBUTION_ANALYSIS': '${ENABLE_CONTRIBUTION_ANALYSIS}',
    'CONTRIBUTION_THRESHOLD': '${CONTRIBUTION_THRESHOLD}',
    'TOP_N_CONTRIBUTIONS': '${TOP_N_CONTRIBUTIONS}',
    'GENERATE_PAIR_CUBES': '${GENERATE_PAIR_CUBES}',
    'MAX_PAIRS_PER_STATE': '${MAX_PAIRS_PER_STATE}',
    'PAIR_CONTRIBUTION_CUTOFF': '${PAIR_CONTRIBUTION_CUTOFF}',
    'GENERATE_GROUND_STATE_DENSITY': '${GENERATE_GROUND_STATE_DENSITY}',
    'GENERATE_ELECTROSTATIC_POTENTIAL': '${GENERATE_ELECTROSTATIC_POTENTIAL}',
    'GENERATE_DEFORMATION_DENSITY': '${GENERATE_DEFORMATION_DENSITY}',
    # Section 7: Grid Settings
    'USE_GRID_RESOLUTION': '${USE_GRID_RESOLUTION}',
    'BOX_MARGIN': '${BOX_MARGIN}',
    'GRID_SPACING': '${GRID_SPACING}',
    # Section 8: Execution Settings
    'NUM_THREADS': '${NUM_THREADS}',
    'VERBOSE_LEVEL': '${VERBOSE_LEVEL}',
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
    'MONOMER_A_ATOMS': '${MONOMER_A_ATOMS}',
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
    echo -e "${YELLOW}Geometry Optimization:${NC}"
    echo "  Optimize:          ${OPTIMISE_GEOMETRY}"
    if [ "$OPTIMISE_GEOMETRY" = True ]; then
        echo "  Max steps:         ${OPT_MAX_STEPS}"
        echo "  Convergence:       ${OPT_CONV_PARAMS}"
    fi
    echo "  Verbose level:     ${VERBOSE_LEVEL}"
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
    # Create output directory first so log file can be saved there
    mkdir -p "${OUTPUT_DIR}"
    
    # Log file path inside output directory
    LOG_FILE_PATH="${OUTPUT_DIR}/${LOG_FILE}"
    TIME_FILE_PATH="${OUTPUT_DIR}/timing_${LOG_FILE%.log}.txt"
    
    echo -e "${GREEN}Starting ${CALC_TYPE} calculation...${NC}"
    echo -e "${BLUE}Command: python3 ${SCRIPT_NAME}${NC}"
    echo -e "${BLUE}Log file: ${LOG_FILE_PATH}${NC}"
    echo -e "${BLUE}Timing file: ${TIME_FILE_PATH}${NC}"
    echo ""
    
    if [ "$RUN_IN_BACKGROUND" = True ]; then
        # Run in background
        if [ "$VERBOSE" = True ]; then
            { time python3 -u "${SCRIPT_NAME}" > "${LOG_FILE_PATH}" 2>&1 ; } 2> "${TIME_FILE_PATH}" &
            PID=$!
            echo -e "${GREEN}✓ Calculation started in background (PID: ${PID})${NC}"
            echo -e "${BLUE}Monitor with: tail -f ${LOG_FILE_PATH}${NC}"
        else
            { time python3 -u "${SCRIPT_NAME}" > "${LOG_FILE_PATH}" 2>&1 ; } 2> "${TIME_FILE_PATH}" &
            echo -e "${GREEN}✓ Calculation started in background (PID: $!)${NC}"
        fi
    else
        # Run in foreground
        if [ "$VERBOSE" = True ]; then
            { time python3 -u "${SCRIPT_NAME}" 2>&1 | tee "${LOG_FILE_PATH}" ; } 2> "${TIME_FILE_PATH}"
            CALC_EXIT_STATUS=${PIPESTATUS[0]}
        else
            { time python3 -u "${SCRIPT_NAME}" > "${LOG_FILE_PATH}" 2>&1 ; } 2> "${TIME_FILE_PATH}"
            CALC_EXIT_STATUS=$?
        fi
        
        # Check exit status
        if [ ${CALC_EXIT_STATUS} -eq 0 ]; then
            echo ""
            echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
            echo -e "${GREEN}✓ Calculation completed successfully!${NC}"
            echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
            echo ""
            echo -e "${BLUE}Output directory: ${OUTPUT_DIR}/${NC}"
            echo -e "${BLUE}Log file: ${LOG_FILE_PATH}${NC}"
            echo -e "${BLUE}Timing file: ${TIME_FILE_PATH}${NC}"
            echo ""
        else
            echo ""
            echo -e "${RED}═══════════════════════════════════════════════════════════════════${NC}"
            echo -e "${RED}✗ Calculation failed!${NC}"
            echo -e "${RED}═══════════════════════════════════════════════════════════════════${NC}"
            echo ""
            echo -e "${YELLOW}Check log file for errors: ${LOG_FILE_PATH}${NC}"
            echo -e "${YELLOW}Timing file: ${TIME_FILE_PATH}${NC}"
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
