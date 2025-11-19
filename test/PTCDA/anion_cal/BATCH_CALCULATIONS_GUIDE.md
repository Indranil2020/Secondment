# 🚀 Batch Calculations Guide - Multiple Charges

## ✅ Problem Solved: No More Overwriting!

### Before (Problem)
```bash
# Run charge 0
CHARGE=0
./run_cal.sh
# Output: output_gpu/

# Run charge +1
CHARGE=1
./run_cal.sh
# Output: output_gpu/  ← OVERWRITES charge 0!

# Run charge -1
CHARGE=-1
./run_cal.sh
# Output: output_gpu/  ← OVERWRITES everything!
```

**Result:** ❌ Only the last calculation's output is kept!

### After (Solution)
```bash
# Run charge 0
CHARGE=0
./run_cal.sh
# Output: output_gpu_charge0/  ✓

# Run charge +1
CHARGE=1
./run_cal.sh
# Output: output_gpu_charge1/  ✓

# Run charge -1
CHARGE=-1
./run_cal.sh
# Output: output_gpu_charge-1/  ✓
```

**Result:** ✅ Each charge has its own directory!

---

## 📁 Output Directory Structure

### Charge-Specific Directories

```
anion_cal/
├── output_gpu_charge0/          # Neutral (charge=0)
│   ├── homo.cube
│   ├── lumo.cube
│   ├── transition_density_state_1.cube
│   ├── nto_state_1.molden
│   └── contribution_tables.txt
│
├── output_gpu_charge1/          # Cation (charge=+1)
│   ├── homo.cube
│   ├── lumo.cube
│   ├── transition_density_state_1.cube
│   ├── nto_state_1_alpha.molden
│   ├── nto_state_1_beta.molden
│   └── contribution_tables.txt
│
├── output_gpu_charge-1/         # Anion (charge=-1)
│   ├── homo.cube
│   ├── lumo.cube
│   ├── transition_density_state_1.cube
│   ├── nto_state_1_alpha.molden
│   ├── nto_state_1_beta.molden
│   └── contribution_tables.txt
│
├── GPU_charge0_20251118_140530.log
├── GPU_charge1_20251118_141245.log
└── GPU_charge-1_20251118_142010.log
```

### CPU vs GPU Directories

```
output_cpu_charge0/              # CPU, neutral
output_cpu_charge1/              # CPU, cation
output_cpu_charge-1/             # CPU, anion

output_gpu_charge0/              # GPU, neutral
output_gpu_charge1/              # GPU, cation
output_gpu_charge-1/             # GPU, anion
```

---

## 🎯 Usage Methods

### Method 1: Manual (One at a Time)

```bash
# Edit run_cal.sh for each charge
vim run_cal.sh
# Change: CHARGE=0
./run_cal.sh

vim run_cal.sh
# Change: CHARGE=1
./run_cal.sh

vim run_cal.sh
# Change: CHARGE=-1
./run_cal.sh
```

### Method 2: Command Line (Quick)

```bash
# Charge 0
sed -i 's/^CHARGE=.*/CHARGE=0/' run_cal.sh && ./run_cal.sh

# Charge +1
sed -i 's/^CHARGE=.*/CHARGE=1/' run_cal.sh && ./run_cal.sh

# Charge -1
sed -i 's/^CHARGE=.*/CHARGE=-1/' run_cal.sh && ./run_cal.sh
```

### Method 3: Batch Script (Automatic) ⭐ **RECOMMENDED**

```bash
# Edit batch configuration
vim run_batch_charges.sh
# Set: CHARGES="0 1 -1"

# Run all charges automatically
./run_batch_charges.sh
```

---

## 🚀 Batch Script (`run_batch_charges.sh`)

### Configuration

```bash
# List of charges to calculate
CHARGES="0 1 -1"              # Space-separated

# Use GPU or CPU?
USE_GPU=true

# Molecule settings (same for all charges)
XYZ_FILE="PTCDA_clean.xyz"
BASIS_SET="6-31g"
XC_FUNCTIONAL="b3lyp"
NUM_EXCITED_STATES=10

# Run in background?
RUN_IN_BACKGROUND=false

# Wait between jobs (seconds)
WAIT_BETWEEN_JOBS=2
```

### Example Output

```
═══════════════════════════════════════════════════════════════════
         BATCH CALCULATION - MULTIPLE CHARGES
═══════════════════════════════════════════════════════════════════

Configuration:
  Molecule:          PTCDA_clean.xyz
  Basis set:         6-31g
  Functional:        b3lyp
  Excited states:    10
  Charges:           0 1 -1
  Mode:              GPU
  Background:        false

═══════════════════════════════════════════════════════════════════

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
JOB 1/3: Charge = 0
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[... calculation runs ...]

✓ Job 1/3 completed successfully

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
JOB 2/3: Charge = 1
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[... calculation runs ...]

✓ Job 2/3 completed successfully

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
JOB 3/3: Charge = -1
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[... calculation runs ...]

✓ Job 3/3 completed successfully

═══════════════════════════════════════════════════════════════════
         BATCH CALCULATION COMPLETE
═══════════════════════════════════════════════════════════════════

Summary:
  Total jobs:        3
  Charges:           0 1 -1

Output directories:
  Charge 0:       output_gpu_charge0/
  Charge 1:       output_gpu_charge1/
  Charge -1:      output_gpu_charge-1/

Log files:
  GPU_charge0_20251118_140530.log
  GPU_charge1_20251118_141245.log
  GPU_charge-1_20251118_142010.log

═══════════════════════════════════════════════════════════════════
```

---

## 📊 Log File Naming

### With Timestamp (Default)
```bash
AUTO_TIMESTAMP=true

# Results:
GPU_charge0_20251118_140530.log
GPU_charge1_20251118_141245.log
GPU_charge-1_20251118_142010.log
```

### Without Timestamp
```bash
AUTO_TIMESTAMP=false

# Results:
GPU_charge0.log
GPU_charge1.log
GPU_charge-1.log
```

---

## 🎯 Common Workflows

### Workflow 1: PTCDA (0, +1, -1)

```bash
# Edit run_batch_charges.sh:
CHARGES="0 1 -1"
XYZ_FILE="PTCDA_clean.xyz"
USE_GPU=true

# Run:
./run_batch_charges.sh

# Results:
# output_gpu_charge0/
# output_gpu_charge1/
# output_gpu_charge-1/
```

### Workflow 2: Test Multiple Charges on CPU First

```bash
# Edit run_batch_charges.sh:
CHARGES="0 1 -1"
USE_GPU=false
NUM_EXCITED_STATES=5
USE_TDA=true

# Run quick test:
./run_batch_charges.sh

# Then full GPU run:
sed -i 's/^USE_GPU=.*/USE_GPU=true/' run_batch_charges.sh
sed -i 's/^NUM_EXCITED_STATES=.*/NUM_EXCITED_STATES=10/' run_batch_charges.sh
sed -i 's/^USE_TDA=.*/USE_TDA=false/' run_batch_charges.sh
./run_batch_charges.sh
```

### Workflow 3: Compare CPU vs GPU (Same Charges)

```bash
# CPU run
CHARGES="0 1 -1"
USE_GPU=false
./run_batch_charges.sh

# GPU run
USE_GPU=true
./run_batch_charges.sh

# Compare results:
for charge in 0 1 -1; do
    echo "Charge ${charge}:"
    diff <(grep "State 1:" CPU_charge${charge}*.log) \
         <(grep "State 1:" GPU_charge${charge}*.log)
done
```

### Workflow 4: Extended Charge Range

```bash
# Edit run_batch_charges.sh:
CHARGES="-2 -1 0 1 2"

# Run:
./run_batch_charges.sh

# Results:
# output_gpu_charge-2/
# output_gpu_charge-1/
# output_gpu_charge0/
# output_gpu_charge1/
# output_gpu_charge2/
```

---

## 📈 Performance Estimates

### PTCDA (3 charges: 0, +1, -1)

| Mode | Time per Charge | Total Time (3 charges) |
|------|-----------------|------------------------|
| **CPU** | ~20-40 min | ~60-120 min (1-2 hours) |
| **GPU** | ~6-10 min | ~18-30 min (< 30 min) |

### Recommendations

**For Testing:**
```bash
USE_GPU=false
USE_TDA=true
NUM_EXCITED_STATES=5
```
**Time:** ~30-60 min for 3 charges

**For Production:**
```bash
USE_GPU=true
USE_TDA=false
NUM_EXCITED_STATES=10
```
**Time:** ~20-30 min for 3 charges

---

## 🔍 Comparing Results

### Quick Comparison Script

```bash
#!/bin/bash
# compare_charges.sh

echo "Comparing results for different charges:"
echo ""

for charge in 0 1 -1; do
    echo "═══════════════════════════════════════"
    echo "Charge: ${charge}"
    echo "═══════════════════════════════════════"
    
    log=$(ls GPU_charge${charge}*.log 2>/dev/null | head -1)
    if [ -f "$log" ]; then
        echo "Ground state energy:"
        grep "Ground state energy:" "$log"
        echo ""
        echo "First 3 excited states:"
        grep "State [1-3]:" "$log" | head -3
    else
        echo "Log file not found"
    fi
    echo ""
done
```

---

## ✅ Advantages

| Feature | Old Method | New Method |
|---------|------------|------------|
| **Overwriting** | ❌ Yes | ✅ No |
| **Organization** | ❌ Poor | ✅ Excellent |
| **Batch Jobs** | ❌ Manual | ✅ Automatic |
| **Log Files** | ❌ Generic | ✅ Charge-specific |
| **Comparison** | ❌ Difficult | ✅ Easy |

---

## 🎉 Summary

### ✅ **Problem Solved!**

**Before:**
- ❌ Outputs overwrite each other
- ❌ Manual tracking needed
- ❌ Difficult to compare

**After:**
- ✅ Each charge has its own directory
- ✅ Automatic organization
- ✅ Easy comparison
- ✅ Batch processing available

### 🚀 **Three Ways to Run Multiple Charges:**

1. **Manual:** Edit `run_cal.sh` for each charge
2. **Command line:** Use `sed` to update and run
3. **Batch:** Use `run_batch_charges.sh` (recommended!)

### 📁 **Output Structure:**

```
output_gpu_charge0/      # Neutral
output_gpu_charge1/      # Cation
output_gpu_charge-1/     # Anion
```

**Everything organized and nothing overwritten!** 🎯
