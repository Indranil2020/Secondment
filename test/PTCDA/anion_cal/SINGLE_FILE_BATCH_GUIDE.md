# 🚀 Single File for Everything - Batch Mode Guide

## ✅ **ONE FILE DOES IT ALL!**

**`run_cal.sh`** now handles both single and batch calculations automatically!

---

## 🎯 How It Works

### Single Charge (Old Way - Still Works)
```bash
# Edit run_cal.sh:
CHARGE=-1

# Run:
./run_cal.sh

# Output:
# output_gpu_charge-1/
# GPU_charge-1_timestamp.log
```

### Multiple Charges (New Batch Mode)
```bash
# Edit run_cal.sh:
CHARGE="0 1 -1"    # ← Just add space-separated values!

# Run:
./run_cal.sh

# Output:
# output_gpu_charge0/
# output_gpu_charge1/
# output_gpu_charge-1/
# GPU_charge0_timestamp.log
# GPU_charge1_timestamp.log
# GPU_charge-1_timestamp.log
```

**That's it!** The script automatically detects batch mode and runs all charges! 🎉

---

## 📊 Examples

### Example 1: PTCDA (Neutral, Cation, Anion)
```bash
# Edit run_cal.sh:
CHARGE="0 1 -1"
XYZ_FILE="PTCDA_clean.xyz"
USE_GPU=true

# Run:
./run_cal.sh
```

**Output:**
```
═══════════════════════════════════════════════════════════════════
         BATCH MODE: 3 CALCULATIONS
═══════════════════════════════════════════════════════════════════

Charges to calculate: 0 1 -1
Mode: GPU

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
  Successful:        3
  Failed:            0

Output directories:
  Charge 0:       output_gpu_charge0/
  Charge 1:       output_gpu_charge1/
  Charge -1:      output_gpu_charge-1/
```

### Example 2: Single Charge (Still Works)
```bash
# Edit run_cal.sh:
CHARGE=-1          # ← Single value (no spaces)

# Run:
./run_cal.sh
```

**Output:** Normal single calculation (no batch mode)

### Example 3: Extended Range
```bash
# Edit run_cal.sh:
CHARGE="-2 -1 0 1 2"

# Run:
./run_cal.sh
```

**Output:** 5 calculations automatically!

---

## 🎛️ Configuration

### In `run_cal.sh`:

```bash
# ============================================================================
# CALCULATION CONFIGURATION
# ============================================================================

# --- Molecule Selection ---
USE_XYZ=true
XYZ_FILE="PTCDA_clean.xyz"
BASIS_SET="6-31g"

# --- Charge and Spin Settings ---
# CHARGE can be:
#   - Single value: CHARGE=-1
#   - Multiple values (batch): CHARGE="0 1 -1"
CHARGE="0 1 -1"              # ← BATCH MODE!

# --- DFT/TDDFT Settings ---
XC_FUNCTIONAL="b3lyp"
NUM_EXCITED_STATES=10
USE_TDA=false

# --- Execution Control ---
USE_GPU=true                 # true=GPU, false=CPU
```

---

## 📁 Output Structure

### Batch Mode: `CHARGE="0 1 -1"`

```
anion_cal/
├── output_gpu_charge0/          # Neutral
│   ├── homo.cube
│   ├── lumo.cube
│   ├── transition_density_state_1.cube
│   └── contribution_tables.txt
│
├── output_gpu_charge1/          # Cation
│   ├── homo.cube
│   ├── lumo.cube
│   ├── nto_state_1_alpha.molden
│   └── nto_state_1_beta.molden
│
├── output_gpu_charge-1/         # Anion
│   ├── homo.cube
│   ├── lumo.cube
│   ├── nto_state_1_alpha.molden
│   └── nto_state_1_beta.molden
│
├── GPU_charge0_20251118_140530.log
├── GPU_charge1_20251118_141245.log
└── GPU_charge-1_20251118_142010.log
```

---

## 🔍 How Detection Works

The script checks if `CHARGE` contains spaces:

```bash
# Single charge (no spaces)
CHARGE=-1          # → Single mode

# Multiple charges (has spaces)
CHARGE="0 1 -1"    # → Batch mode
CHARGE="0 1"       # → Batch mode
CHARGE="-1 0 1 2"  # → Batch mode
```

---

## ✅ Advantages

### Before (Separate Batch Script)
```
run_cal.sh              # Single calculations
run_batch_charges.sh    # Batch calculations
```
**Problem:** Two files to maintain ❌

### After (Integrated)
```
run_cal.sh              # Single AND batch!
```
**Solution:** One file for everything! ✅

---

## 🎯 Usage Comparison

| Mode | Configuration | Command |
|------|---------------|---------|
| **Single** | `CHARGE=-1` | `./run_cal.sh` |
| **Batch** | `CHARGE="0 1 -1"` | `./run_cal.sh` |

**Same command, different modes!** 🚀

---

## 📊 Quick Reference

### Single Charge
```bash
CHARGE=-1              # No quotes needed
```

### Batch Charges
```bash
CHARGE="0 1 -1"        # Quotes + spaces
CHARGE="0 1"           # Any number of charges
CHARGE="-2 -1 0 1 2"   # Extended range
```

### Switch CPU/GPU
```bash
USE_GPU=true           # GPU (faster)
USE_GPU=false          # CPU (slower)
```

---

## 🎉 Summary

### ✅ **ONE FILE FOR EVERYTHING!**

**Single Charge:**
```bash
CHARGE=-1
./run_cal.sh
```

**Multiple Charges:**
```bash
CHARGE="0 1 -1"
./run_cal.sh
```

**That's it!** No separate batch script needed! 🎯

### Key Points:
- ✅ Automatic batch detection
- ✅ Charge-specific directories
- ✅ Timestamped log files
- ✅ Progress tracking
- ✅ Summary at the end
- ✅ One file to maintain

**Simpler, cleaner, better!** 🚀
