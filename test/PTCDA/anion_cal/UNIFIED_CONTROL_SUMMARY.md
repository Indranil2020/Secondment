# 🚀 Unified Calculation Control - Summary

## ✅ **NEW: `run_cal.sh` - One Script for Everything!**

**Date:** November 18, 2025  
**Status:** ✅ READY TO USE

---

## 🎯 What's New

### Before (Old Method)
```bash
# Edit Python script manually
vim tdm_calc_accurate_GPU.py
# Change CHARGE = -1
# Change XYZ_FILE = "molecule.xyz"
# etc...

# Run
./run_gpu.sh

# For CPU, edit different script
vim tdm_calc_accurate_cpu.py
# Change same settings again...
python3 tdm_calc_accurate_cpu.py
```

**Problems:**
- ❌ Edit Python scripts directly
- ❌ Different process for CPU/GPU
- ❌ Settings scattered in code
- ❌ Error-prone manual editing
- ❌ No validation

### After (New Method)
```bash
# Edit ONE bash script
vim run_cal.sh
# Change CHARGE=-1
# Change XYZ_FILE="molecule.xyz"
# Change USE_GPU=true  # or false for CPU

# Run (same command for CPU or GPU!)
./run_cal.sh
```

**Benefits:**
- ✅ Edit ONE control script
- ✅ Same process for CPU/GPU
- ✅ All settings in one place
- ✅ Automatic validation
- ✅ Professional output

---

## 📊 Quick Comparison

| Feature | Old Method | New Method (`run_cal.sh`) |
|---------|------------|---------------------------|
| **Configuration** | Edit Python | Edit bash |
| **CPU/GPU Switch** | Different scripts | `USE_GPU=true/false` |
| **Settings Location** | Scattered in code | Top of script |
| **Validation** | None | Automatic |
| **Error Messages** | Basic | Colored, detailed |
| **Log Files** | Manual | Auto-timestamped |
| **Learning Curve** | Need Python | Just bash variables |
| **Batch Jobs** | Complex | Simple |

---

## 🎛️ Configuration Example

### All Settings in One Place
```bash
# ============================================================================
# CALCULATION CONFIGURATION (Top of run_cal.sh)
# ============================================================================

# Molecule
USE_XYZ=true
XYZ_FILE="PTCDA_clean.xyz"
BASIS_SET="6-31g"
CHARGE=-1
SPIN=None

# DFT/TDDFT
XC_FUNCTIONAL="b3lyp"
NUM_EXCITED_STATES=10
USE_TDA=false

# Output
STATES_TO_OUTPUT="0 1 2"
GENERATE_TRANSITION_DENSITY=true
GENERATE_EXCITED_DENSITY=true
GENERATE_DENSITY_DIFFERENCE=true
GENERATE_HOMO_LUMO=true

# Analysis
ENABLE_NTO_ANALYSIS=true
NTO_STATES="0 1 2"
ENABLE_CONTRIBUTION_ANALYSIS=true
CONTRIBUTION_STATES="0 1 2"

# Execution
USE_GPU=true                    # ← SINGLE SWITCH FOR CPU/GPU!
LOG_FILE="calculation.log"
RUN_IN_BACKGROUND=false
```

---

## 🚀 Usage Examples

### Example 1: Run PTCDA Anion on GPU
```bash
# Edit run_cal.sh:
USE_GPU=true
XYZ_FILE="PTCDA_clean.xyz"
CHARGE=-1

# Run:
./run_cal.sh
```

### Example 2: Switch to CPU (Same Settings)
```bash
# Edit run_cal.sh:
USE_GPU=false    # ← ONLY CHANGE THIS!

# Run:
./run_cal.sh
```

### Example 3: Test Different Charges
```bash
# Neutral
CHARGE=0
./run_cal.sh

# Cation
CHARGE=1
./run_cal.sh

# Anion
CHARGE=-1
./run_cal.sh
```

---

## 🎨 Professional Output

### Configuration Display
```
═══════════════════════════════════════════════════════════════════
         TDDFT CALCULATION CONFIGURATION
═══════════════════════════════════════════════════════════════════

Execution:
  Mode:              GPU
  Script:            tdm_calc_accurate_GPU.py
  Output directory:  output_gpu
  Log file:          GPU_20251118_135900.log

Molecule:
  Use XYZ:           true
  XYZ file:          PTCDA_clean.xyz
  Basis set:         6-31g
  Charge:            -1
  Spin:              None

DFT/TDDFT:
  Functional:        b3lyp
  Excited states:    10
  Use TDA:           false

═══════════════════════════════════════════════════════════════════

Updating configuration in tdm_calc_accurate_GPU.py...
✓ Configuration updated successfully

Starting GPU calculation...
Command: python3 tdm_calc_accurate_GPU.py
Log file: GPU_20251118_135900.log

[... calculation runs ...]

═══════════════════════════════════════════════════════════════════
✓ Calculation completed successfully!
═══════════════════════════════════════════════════════════════════

Output directory: output_gpu/
Log file: GPU_20251118_135900.log
```

---

## 📁 File Structure

### Before
```
anion_cal/
├── tdm_calc_accurate_cpu.py      # Edit this for CPU
├── tdm_calc_accurate_GPU.py      # Edit this for GPU
├── run_gpu.sh                    # Run GPU
└── (no unified control)
```

### After
```
anion_cal/
├── tdm_calc_accurate_cpu.py      # Don't edit directly
├── tdm_calc_accurate_GPU.py      # Don't edit directly
├── run_cal.sh                    # ← EDIT THIS! (unified control)
├── run_gpu.sh                    # (legacy, can keep)
├── RUN_CAL_GUIDE.md              # User guide
└── UNIFIED_CONTROL_SUMMARY.md    # This file
```

---

## 🔧 How It Works

### 1. You Edit `run_cal.sh`
```bash
CHARGE=-1
USE_GPU=true
```

### 2. Script Updates Python File
```bash
# Automatically updates tdm_calc_accurate_GPU.py:
CHARGE = -1
```

### 3. Script Runs Calculation
```bash
python3 tdm_calc_accurate_GPU.py > GPU_20251118_135900.log 2>&1
```

### 4. Reports Results
```bash
✓ Calculation completed successfully!
Output directory: output_gpu/
```

---

## 🎯 Common Workflows

### Workflow 1: Batch Testing
```bash
# Test multiple charges
for charge in 0 1 -1; do
    sed -i "s/^CHARGE=.*/CHARGE=$charge/" run_cal.sh
    ./run_cal.sh
done
```

### Workflow 2: CPU/GPU Comparison
```bash
# GPU run
sed -i 's/^USE_GPU=.*/USE_GPU=true/' run_cal.sh
./run_cal.sh

# CPU run (same settings)
sed -i 's/^USE_GPU=.*/USE_GPU=false/' run_cal.sh
./run_cal.sh

# Compare
diff <(grep "State 1:" GPU_*.log) <(grep "State 1:" CPU_*.log)
```

### Workflow 3: Progressive Analysis
```bash
# Quick test
sed -i 's/^USE_TDA=.*/USE_TDA=true/' run_cal.sh
sed -i 's/^NUM_EXCITED_STATES=.*/NUM_EXCITED_STATES=5/' run_cal.sh
./run_cal.sh

# Full analysis
sed -i 's/^USE_TDA=.*/USE_TDA=false/' run_cal.sh
sed -i 's/^NUM_EXCITED_STATES=.*/NUM_EXCITED_STATES=10/' run_cal.sh
./run_cal.sh
```

---

## ✅ Key Features

### 1. **Centralized Configuration**
All settings in one place at the top of `run_cal.sh`

### 2. **CPU/GPU Switch**
```bash
USE_GPU=true   # GPU
USE_GPU=false  # CPU
```

### 3. **Automatic Updates**
Script updates Python files automatically - no manual editing!

### 4. **Validation**
Checks for file existence, valid settings, etc.

### 5. **Professional Output**
- Colored output (green=success, red=error, blue=info)
- Formatted configuration display
- Progress tracking

### 6. **Timestamped Logs**
```bash
AUTO_TIMESTAMP=true
# Creates: GPU_20251118_135900.log
```

### 7. **Background Execution**
```bash
RUN_IN_BACKGROUND=true
# Run in background, monitor with: tail -f GPU_*.log
```

---

## 📊 Performance

### No Overhead
- Configuration update: < 0.1 seconds
- Same calculation speed as before
- Just better control and logging

---

## 🎉 Summary

### ✅ **One Script to Control Everything!**

**Before:**
- Edit Python scripts manually
- Different process for CPU/GPU
- Error-prone
- No validation

**After:**
- Edit ONE bash script
- Same process for CPU/GPU
- Automatic validation
- Professional output

### 🚀 **Ready for Production!**

```bash
# Just edit and run!
vim run_cal.sh
./run_cal.sh
```

**That's it!** 🎯

---

## 📚 Documentation

1. ✅ `run_cal.sh` - The unified control script
2. ✅ `RUN_CAL_GUIDE.md` - Detailed user guide
3. ✅ `UNIFIED_CONTROL_SUMMARY.md` - This summary

**Everything you need to run calculations efficiently!** 🚀
