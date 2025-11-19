# 🚀 Unified Calculation Control Script - User Guide

## 📋 Overview

**`run_cal.sh`** - Single unified script to control all TDDFT calculations (CPU or GPU)

**Key Features:**
- ✅ Centralized configuration (all settings in one place)
- ✅ CPU/GPU selection with single switch
- ✅ Automatic script configuration
- ✅ Colored output and progress tracking
- ✅ Timestamped log files
- ✅ Background/foreground execution

---

## 🎯 Quick Start

### 1. Basic Usage (GPU, Default Settings)
```bash
./run_cal.sh
```

### 2. Switch to CPU
```bash
# Edit run_cal.sh, change:
USE_GPU=false

# Then run:
./run_cal.sh
```

### 3. Change Molecule
```bash
# Edit run_cal.sh, change:
XYZ_FILE="my_molecule.xyz"
CHARGE=0
SPIN=None

# Then run:
./run_cal.sh
```

---

## ⚙️ Configuration Sections

### 1. **Molecule Selection**
```bash
USE_XYZ=true                    # Use XYZ file or test molecule
XYZ_FILE="PTCDA_clean.xyz"      # Path to XYZ file
BASIS_SET="6-31g"               # Basis set
```

### 2. **Charge and Spin**
```bash
CHARGE=-1                       # 0=neutral, +1=cation, -1=anion
SPIN=None                       # None=auto, 1=singlet, 2=doublet
```

### 3. **DFT/TDDFT Settings**
```bash
XC_FUNCTIONAL="b3lyp"           # b3lyp, pbe0, cam-b3lyp, etc.
NUM_EXCITED_STATES=10           # Number of states
USE_TDA=false                   # true=TDA (fast), false=TDDFT (accurate)
```

### 4. **Output Selection**
```bash
STATES_TO_OUTPUT="0 1 2"        # Which states for cube files
GENERATE_TRANSITION_DENSITY=true
GENERATE_EXCITED_DENSITY=true
GENERATE_DENSITY_DIFFERENCE=true
GENERATE_HOMO_LUMO=true
```

### 5. **Grid Settings**
```bash
USE_GRID_RESOLUTION=false       # Fixed resolution or spacing
GRID_RESOLUTION_X=80            # Grid points (if fixed)
GRID_RESOLUTION_Y=80
GRID_RESOLUTION_Z=80
BOX_MARGIN=4.0                  # Margin in Angstrom
GRID_SPACING=0.2                # Spacing in Angstrom
```

### 6. **NTO Analysis**
```bash
ENABLE_NTO_ANALYSIS=true        # Enable NTO analysis
NTO_STATES="0 1 2"              # Which states
```

### 7. **Contribution Analysis**
```bash
ENABLE_CONTRIBUTION_ANALYSIS=true
CONTRIBUTION_STATES="0 1 2"
TOP_N_CONTRIBUTIONS=10
GENERATE_PAIR_CUBES=true
MAX_PAIRS_PER_STATE=3
PAIR_CONTRIBUTION_CUTOFF=0.05
```

### 8. **Execution Control**
```bash
USE_GPU=true                    # true=GPU, false=CPU
LOG_FILE="calculation.log"      # Log file name
AUTO_TIMESTAMP=true             # Add timestamp to log
RUN_IN_BACKGROUND=false         # Background execution
VERBOSE=true                    # Show progress
```

---

## 📊 Usage Examples

### Example 1: PTCDA Anion on GPU
```bash
# Edit run_cal.sh:
USE_GPU=true
XYZ_FILE="PTCDA_clean.xyz"
CHARGE=-1
SPIN=None
XC_FUNCTIONAL="b3lyp"
NUM_EXCITED_STATES=10
STATES_TO_OUTPUT="0 1 2"

# Run:
./run_cal.sh
```

### Example 2: H2O Neutral on CPU
```bash
# Edit run_cal.sh:
USE_GPU=false
XYZ_FILE="H2O.xyz"
CHARGE=0
SPIN=None
NUM_EXCITED_STATES=5
STATES_TO_OUTPUT="0"

# Run:
./run_cal.sh
```

### Example 3: Quick Test (TDA, Few States)
```bash
# Edit run_cal.sh:
USE_GPU=true
USE_TDA=true                    # Faster
NUM_EXCITED_STATES=5            # Fewer states
STATES_TO_OUTPUT="0"            # Only first state
GENERATE_PAIR_CUBES=false       # Skip pair cubes

# Run:
./run_cal.sh
```

### Example 4: Detailed Analysis (All Features)
```bash
# Edit run_cal.sh:
USE_GPU=true
NUM_EXCITED_STATES=10
STATES_TO_OUTPUT="0 1 2 3 4"
ENABLE_NTO_ANALYSIS=true
NTO_STATES="0 1 2 3 4"
ENABLE_CONTRIBUTION_ANALYSIS=true
CONTRIBUTION_STATES="0 1 2 3 4"
GENERATE_PAIR_CUBES=true
MAX_PAIRS_PER_STATE=5

# Run:
./run_cal.sh
```

### Example 5: Background Execution
```bash
# Edit run_cal.sh:
RUN_IN_BACKGROUND=true
VERBOSE=true

# Run:
./run_cal.sh

# Monitor:
tail -f GPU_20251118_135900.log
```

---

## 🎨 Output Display

### Configuration Summary
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

Output:
  States for cubes:  [0 1 2]
  Transition dens:   true
  Excited density:   true
  Density diff:      true
  HOMO/LUMO:         true

Analysis:
  NTO analysis:      true
    NTO states:      [0 1 2]
  Contribution:      true
    States:          [0 1 2]
    Top N pairs:     10
    Pair cubes:      true

═══════════════════════════════════════════════════════════════════
```

---

## 📁 Output Files

### Log Files
```
GPU_20251118_135900.log         # Timestamped log (if AUTO_TIMESTAMP=true)
calculation.log                 # Fixed name (if AUTO_TIMESTAMP=false)
```

### Output Directories
```
output/                         # CPU output
output_gpu/                     # GPU output
```

---

## 🔧 How It Works

### 1. **Configuration Update**
The script automatically updates the Python script with your settings:
```bash
# You set in run_cal.sh:
CHARGE=-1

# Script updates in tdm_calc_accurate_GPU.py:
CHARGE = -1
```

### 2. **List Conversion**
Space-separated lists are converted to Python lists:
```bash
# You set:
STATES_TO_OUTPUT="0 1 2"

# Becomes:
STATES_TO_OUTPUT = [0, 1, 2]
```

### 3. **Execution**
- Updates configuration
- Runs Python script
- Captures output to log file
- Reports success/failure

---

## 🎯 Common Workflows

### Workflow 1: Test Different Charges
```bash
# Test neutral
CHARGE=0
./run_cal.sh

# Test cation
CHARGE=1
./run_cal.sh

# Test anion
CHARGE=-1
./run_cal.sh
```

### Workflow 2: CPU vs GPU Comparison
```bash
# Run on GPU
USE_GPU=true
./run_cal.sh

# Run on CPU
USE_GPU=false
./run_cal.sh

# Compare results
diff <(grep "State 1:" GPU_*.log) <(grep "State 1:" CPU_*.log)
```

### Workflow 3: Progressive Analysis
```bash
# Step 1: Quick test (TDA, 5 states)
USE_TDA=true
NUM_EXCITED_STATES=5
STATES_TO_OUTPUT="0"
./run_cal.sh

# Step 2: Full TDDFT (10 states)
USE_TDA=false
NUM_EXCITED_STATES=10
STATES_TO_OUTPUT="0 1 2"
./run_cal.sh

# Step 3: Detailed analysis (all features)
ENABLE_NTO_ANALYSIS=true
ENABLE_CONTRIBUTION_ANALYSIS=true
./run_cal.sh
```

---

## 🚨 Troubleshooting

### Issue 1: Script Not Found
```bash
# Error: tdm_calc_accurate_GPU.py not found
# Solution: Check you're in the correct directory
cd /home/indranil/Documents/Secondment/test/PTCDA/anion_cal
./run_cal.sh
```

### Issue 2: XYZ File Not Found
```bash
# Error: Cannot open XYZ file
# Solution: Check XYZ_FILE path is correct
XYZ_FILE="PTCDA_clean.xyz"  # Relative path
# or
XYZ_FILE="/full/path/to/PTCDA_clean.xyz"  # Absolute path
```

### Issue 3: GPU Not Available
```bash
# Error: CUDA not found
# Solution: Switch to CPU
USE_GPU=false
./run_cal.sh
```

### Issue 4: Configuration Not Updated
```bash
# If changes don't take effect:
# 1. Check you saved run_cal.sh
# 2. Check Python script was updated
grep "CHARGE =" tdm_calc_accurate_GPU.py
```

---

## 📊 Performance Tips

### For Large Molecules (PTCDA)
```bash
USE_GPU=true                    # Use GPU for speed
USE_TDA=false                   # Full TDDFT for accuracy
NUM_EXCITED_STATES=10           # Reasonable number
STATES_TO_OUTPUT="0 1 2"        # Limit cube files
GENERATE_PAIR_CUBES=true        # Useful analysis
MAX_PAIRS_PER_STATE=3           # Limit to top 3
```

### For Quick Tests
```bash
USE_GPU=true                    # GPU still faster
USE_TDA=true                    # 2× faster
NUM_EXCITED_STATES=5            # Fewer states
STATES_TO_OUTPUT="0"            # Only first state
GENERATE_PAIR_CUBES=false       # Skip to save time
```

### For Detailed Analysis
```bash
USE_GPU=true                    # Speed
NUM_EXCITED_STATES=10           # More states
STATES_TO_OUTPUT="0 1 2 3 4"    # More cubes
ENABLE_NTO_ANALYSIS=true        # Full NTO
ENABLE_CONTRIBUTION_ANALYSIS=true
MAX_PAIRS_PER_STATE=5           # More pairs
```

---

## ✅ Advantages Over Old Method

| Feature | Old (run_gpu.sh) | New (run_cal.sh) |
|---------|------------------|------------------|
| **Configuration** | Edit Python script | Edit bash script |
| **CPU/GPU Switch** | Different scripts | Single switch |
| **Settings** | Scattered | Centralized |
| **Validation** | Manual | Automatic |
| **Log Files** | Manual naming | Auto timestamp |
| **Progress** | Basic | Colored output |
| **Error Handling** | Limited | Comprehensive |

---

## 🎉 Summary

**One Script to Rule Them All!** 🚀

✅ **Centralized control** - All settings in one place  
✅ **Easy switching** - CPU/GPU with one variable  
✅ **Automatic updates** - No manual Python editing  
✅ **Professional output** - Colored, formatted display  
✅ **Flexible execution** - Foreground/background  
✅ **Production ready** - Error handling, logging  

**Just edit `run_cal.sh` and run!** 🎯
