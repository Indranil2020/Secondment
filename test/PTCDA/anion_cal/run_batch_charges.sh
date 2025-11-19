#!/bin/bash
################################################################################
# BATCH CALCULATION SCRIPT - Multiple Charges
################################################################################
# This script runs TDDFT calculations for multiple charge states automatically.
# Each charge gets its own output directory and log file.
#
# Example: PTCDA with charges 0, +1, -1
################################################################################

# ============================================================================
# BATCH CONFIGURATION
# ============================================================================

# List of charges to calculate (space-separated)
CHARGES="0 1 -1"

# Use GPU or CPU?
USE_GPU=true

# Molecule settings (same for all charges)
XYZ_FILE="PTCDA_clean.xyz"
BASIS_SET="6-31g"
XC_FUNCTIONAL="b3lyp"
NUM_EXCITED_STATES=10

# Run in background? (recommended for batch jobs)
RUN_IN_BACKGROUND=false

# Wait between jobs (seconds) - useful if running in foreground
WAIT_BETWEEN_JOBS=2

# ============================================================================
# EXECUTION
# ============================================================================

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}         BATCH CALCULATION - MULTIPLE CHARGES${NC}"
echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
echo ""
echo -e "${YELLOW}Configuration:${NC}"
echo "  Molecule:          ${XYZ_FILE}"
echo "  Basis set:         ${BASIS_SET}"
echo "  Functional:        ${XC_FUNCTIONAL}"
echo "  Excited states:    ${NUM_EXCITED_STATES}"
echo "  Charges:           ${CHARGES}"
echo "  Mode:              $([ "$USE_GPU" = true ] && echo "GPU" || echo "CPU")"
echo "  Background:        ${RUN_IN_BACKGROUND}"
echo ""
echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
echo ""

# Counter
TOTAL_JOBS=$(echo $CHARGES | wc -w)
CURRENT_JOB=0

# Loop through charges
for CHARGE in $CHARGES; do
    CURRENT_JOB=$((CURRENT_JOB + 1))
    
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}JOB ${CURRENT_JOB}/${TOTAL_JOBS}: Charge = ${CHARGE}${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    
    # Update run_cal.sh with current charge
    sed -i "s/^CHARGE=.*/CHARGE=${CHARGE}/" run_cal.sh
    sed -i "s/^USE_GPU=.*/USE_GPU=${USE_GPU}/" run_cal.sh
    sed -i "s/^XYZ_FILE=.*/XYZ_FILE=\"${XYZ_FILE}\"/" run_cal.sh
    sed -i "s/^BASIS_SET=.*/BASIS_SET=\"${BASIS_SET}\"/" run_cal.sh
    sed -i "s/^XC_FUNCTIONAL=.*/XC_FUNCTIONAL=\"${XC_FUNCTIONAL}\"/" run_cal.sh
    sed -i "s/^NUM_EXCITED_STATES=.*/NUM_EXCITED_STATES=${NUM_EXCITED_STATES}/" run_cal.sh
    sed -i "s/^RUN_IN_BACKGROUND=.*/RUN_IN_BACKGROUND=${RUN_IN_BACKGROUND}/" run_cal.sh
    
    # Run calculation
    ./run_cal.sh
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Job ${CURRENT_JOB}/${TOTAL_JOBS} completed successfully${NC}"
    else
        echo -e "${RED}✗ Job ${CURRENT_JOB}/${TOTAL_JOBS} failed${NC}"
    fi
    
    # Wait between jobs (if not last job and not in background)
    if [ $CURRENT_JOB -lt $TOTAL_JOBS ] && [ "$RUN_IN_BACKGROUND" = false ]; then
        echo ""
        echo -e "${YELLOW}Waiting ${WAIT_BETWEEN_JOBS} seconds before next job...${NC}"
        sleep $WAIT_BETWEEN_JOBS
    fi
done

echo ""
echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}         BATCH CALCULATION COMPLETE${NC}"
echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
echo ""
echo -e "${YELLOW}Summary:${NC}"
echo "  Total jobs:        ${TOTAL_JOBS}"
echo "  Charges:           ${CHARGES}"
echo ""
echo -e "${YELLOW}Output directories:${NC}"
for CHARGE in $CHARGES; do
    if [ "$CHARGE" -ge 0 ]; then
        CHARGE_SUFFIX="charge${CHARGE}"
    else
        CHARGE_SUFFIX="charge${CHARGE}"
    fi
    
    if [ "$USE_GPU" = true ]; then
        OUTPUT_DIR="output_gpu_${CHARGE_SUFFIX}"
    else
        OUTPUT_DIR="output_cpu_${CHARGE_SUFFIX}"
    fi
    
    echo "  Charge ${CHARGE}:       ${OUTPUT_DIR}/"
done
echo ""
echo -e "${YELLOW}Log files:${NC}"
ls -1 *charge*.log 2>/dev/null | tail -${TOTAL_JOBS} | while read log; do
    echo "  ${log}"
done
echo ""
echo -e "${GREEN}═══════════════════════════════════════════════════════════════════${NC}"
