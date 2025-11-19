# Configuration Guide for tdm_calc_accurate.py

## Overview

The script now includes a comprehensive configuration section at the top of the file (lines 22-66) that allows you to control all aspects of the calculation without modifying the core code.

## Configuration Options

### 1. Parallel Calculation Settings

```python
ENABLE_PARALLEL = True  # Enable/disable parallel computation
NUM_THREADS = 4         # Number of CPU threads (set to 0 for auto-detect)
```

**Options:**
- `ENABLE_PARALLEL = True`: Use multiple CPU cores
- `ENABLE_PARALLEL = False`: Use single core (useful for debugging)
- `NUM_THREADS = 4`: Use 4 threads
- `NUM_THREADS = 0`: Auto-detect and use all available cores

**Recommendation:** 
- For PTCDA: Use 4-8 threads for optimal performance
- For small molecules: 2-4 threads is sufficient
- Set to 0 for automatic detection

---

### 2. Molecule Selection

```python
USE_XYZ = True                    # True: use XYZ file, False: use H2O
XYZ_FILE = 'PTCDA.xyz'           # Path to XYZ file
BASIS_SET = '6-31g'              # Basis set
```

**Options:**
- `USE_XYZ = False`: Use built-in H2O test molecule
- `USE_XYZ = True`: Load molecule from XYZ file
- `XYZ_FILE`: Path to your XYZ file (relative or absolute)
- `BASIS_SET`: Any PySCF-supported basis set
  - `'6-31g'`: Fast, good for testing
  - `'6-31g*'` or `'6-31g(d)'`: Better for excited states
  - `'def2-SVP'`: Good balance
  - `'def2-TZVP'`: High quality (slower)

---

### 3. TDDFT Settings

```python
NUM_EXCITED_STATES = 10  # Total number of excited states to calculate
```

**Options:**
- Set to any integer (typically 5-20)
- More states = longer calculation time
- You can calculate 10 states but only output cube files for selected ones

**Recommendation:**
- Start with 10 states
- Increase if you need higher excited states

---

### 4. Output Selection

```python
STATES_TO_OUTPUT = [0, 1, 2]  # Which states to generate cube files for (0-indexed)
```

**Important:** States are 0-indexed!
- `[0, 1, 2]` = States 1, 2, 3
- `[0, 4, 9]` = States 1, 5, 10
- `range(5)` = First 5 states (0-4)
- `[0]` = Only first excited state

**Examples:**
```python
# Only first excited state
STATES_TO_OUTPUT = [0]

# First three states
STATES_TO_OUTPUT = [0, 1, 2]

# First five states
STATES_TO_OUTPUT = list(range(5))

# Specific states (1st, 3rd, 5th)
STATES_TO_OUTPUT = [0, 2, 4]

# All calculated states
STATES_TO_OUTPUT = list(range(NUM_EXCITED_STATES))
```

---

### 5. Cube File Generation Options

```python
GENERATE_TRANSITION_DENSITY = True   # Transition density matrix
GENERATE_EXCITED_DENSITY = True      # Excited state density
GENERATE_DENSITY_DIFFERENCE = True   # Density difference (excited - ground)
GENERATE_HOMO_LUMO = True            # HOMO and LUMO orbitals
```

**Turn on/off specific file types:**

| Option | Description | Use Case |
|--------|-------------|----------|
| `GENERATE_TRANSITION_DENSITY` | Transition density matrix | Transition dipole analysis |
| `GENERATE_EXCITED_DENSITY` | Total excited state density | Compare with ground state |
| `GENERATE_DENSITY_DIFFERENCE` | Density redistribution | **Best for visualization** |
| `GENERATE_HOMO_LUMO` | HOMO/LUMO orbitals | Verify transition character |

**Common configurations:**

```python
# Minimal output (only density difference)
GENERATE_TRANSITION_DENSITY = False
GENERATE_EXCITED_DENSITY = False
GENERATE_DENSITY_DIFFERENCE = True
GENERATE_HOMO_LUMO = False

# Full output (all files)
GENERATE_TRANSITION_DENSITY = True
GENERATE_EXCITED_DENSITY = True
GENERATE_DENSITY_DIFFERENCE = True
GENERATE_HOMO_LUMO = True

# Only for verification
GENERATE_TRANSITION_DENSITY = True
GENERATE_EXCITED_DENSITY = False
GENERATE_DENSITY_DIFFERENCE = False
GENERATE_HOMO_LUMO = True
```

---

### 6. Grid Settings

#### Option A: Fixed Grid Resolution (Recommended)

```python
USE_GRID_RESOLUTION = True
GRID_RESOLUTION = [80, 80, 80]  # [nx, ny, nz] grid points
```

**Common resolutions:**
- `[60, 60, 60]`: Fast, lower quality (216,000 points)
- `[80, 80, 80]`: Good balance (512,000 points) ⭐ **Recommended**
- `[100, 100, 100]`: High quality (1,000,000 points)
- `[120, 120, 120]`: Very high quality (1,728,000 points)

**For non-cubic molecules (like PTCDA):**
```python
# PTCDA is flat, so can use different resolutions
GRID_RESOLUTION = [100, 100, 60]  # More points in xy plane, fewer in z
```

#### Option B: Box Dimensions

```python
USE_GRID_RESOLUTION = False
BOX_MARGIN = 4.0      # Margin around molecule in Angstrom
GRID_SPACING = 0.2    # Grid spacing in Angstrom
```

**Settings:**
- `BOX_MARGIN`: Extra space around molecule (3-5 Å typical)
- `GRID_SPACING`: Distance between grid points (0.1-0.3 Å typical)
  - `0.1 Å`: Very fine (slow)
  - `0.2 Å`: Good quality ⭐ **Recommended**
  - `0.3 Å`: Fast, lower quality

**When to use:**
- Use Option A (fixed resolution) for consistent quality
- Use Option B (box dimensions) for automatic sizing

---

### 7. NTO Analysis

```python
ENABLE_NTO_ANALYSIS = True   # Generate NTO molden files
NTO_STATES = [0, 1, 2]       # Which states to perform NTO analysis on
```

**Options:**
- `ENABLE_NTO_ANALYSIS = True`: Generate NTO files
- `ENABLE_NTO_ANALYSIS = False`: Skip NTO analysis (faster)
- `NTO_STATES`: List of states (0-indexed) for NTO analysis

**Note:** NTO analysis is relatively fast, so usually keep it enabled.

---

### 8. Output Directory

```python
OUTPUT_DIR = 'output'  # Directory for output files
```

**Options:**
- `'output'`: Default directory
- `'results'`: Alternative name
- `'ptcda_results'`: Descriptive name
- `'/path/to/results'`: Absolute path

The directory will be created automatically if it doesn't exist.

---

## Example Configurations

### Configuration 1: Quick Test (H2O)

```python
ENABLE_PARALLEL = True
NUM_THREADS = 2
USE_XYZ = False  # Use H2O
NUM_EXCITED_STATES = 5
STATES_TO_OUTPUT = [0, 1]
GENERATE_TRANSITION_DENSITY = True
GENERATE_EXCITED_DENSITY = False
GENERATE_DENSITY_DIFFERENCE = True
GENERATE_HOMO_LUMO = True
USE_GRID_RESOLUTION = True
GRID_RESOLUTION = [60, 60, 60]
ENABLE_NTO_ANALYSIS = True
NTO_STATES = [0, 1]
OUTPUT_DIR = 'test_output'
```

### Configuration 2: Full PTCDA Calculation

```python
ENABLE_PARALLEL = True
NUM_THREADS = 8
USE_XYZ = True
XYZ_FILE = 'PTCDA.xyz'
BASIS_SET = '6-31g'
NUM_EXCITED_STATES = 10
STATES_TO_OUTPUT = [0, 1, 2, 3, 4]  # First 5 states
GENERATE_TRANSITION_DENSITY = True
GENERATE_EXCITED_DENSITY = True
GENERATE_DENSITY_DIFFERENCE = True
GENERATE_HOMO_LUMO = True
USE_GRID_RESOLUTION = True
GRID_RESOLUTION = [100, 100, 60]  # Higher resolution in xy plane
ENABLE_NTO_ANALYSIS = True
NTO_STATES = [0, 1, 2]
OUTPUT_DIR = 'ptcda_results'
```

### Configuration 3: Minimal Output (Only Density Difference)

```python
ENABLE_PARALLEL = True
NUM_THREADS = 4
USE_XYZ = True
XYZ_FILE = 'PTCDA.xyz'
BASIS_SET = '6-31g'
NUM_EXCITED_STATES = 10
STATES_TO_OUTPUT = [0, 1, 2]
GENERATE_TRANSITION_DENSITY = False
GENERATE_EXCITED_DENSITY = False
GENERATE_DENSITY_DIFFERENCE = True  # Only this
GENERATE_HOMO_LUMO = False
USE_GRID_RESOLUTION = True
GRID_RESOLUTION = [80, 80, 80]
ENABLE_NTO_ANALYSIS = False
OUTPUT_DIR = 'output'
```

### Configuration 4: High Quality for Publication

```python
ENABLE_PARALLEL = True
NUM_THREADS = 0  # Auto-detect
USE_XYZ = True
XYZ_FILE = 'PTCDA.xyz'
BASIS_SET = '6-31g*'  # Better basis set
NUM_EXCITED_STATES = 15
STATES_TO_OUTPUT = [0, 1, 2, 3, 4, 5]
GENERATE_TRANSITION_DENSITY = True
GENERATE_EXCITED_DENSITY = True
GENERATE_DENSITY_DIFFERENCE = True
GENERATE_HOMO_LUMO = True
USE_GRID_RESOLUTION = True
GRID_RESOLUTION = [120, 120, 80]  # High resolution
ENABLE_NTO_ANALYSIS = True
NTO_STATES = [0, 1, 2, 3, 4, 5]
OUTPUT_DIR = 'publication_results'
```

---

## Performance Tips

### Calculation Time Estimates (PTCDA, 6-31g basis)

| Configuration | Threads | Grid | Time |
|---------------|---------|------|------|
| Quick test | 2 | 60³ | ~10 min |
| Standard | 4 | 80³ | ~20 min |
| High quality | 8 | 100³ | ~40 min |
| Publication | 8 | 120³ | ~60 min |

### Optimization Tips

1. **Start small:** Test with H2O first
2. **Use parallel:** Always enable parallel computation
3. **Grid resolution:** 80³ is usually sufficient
4. **Selective output:** Only generate files you need
5. **Basis set:** Start with 6-31g, upgrade to 6-31g* if needed

### Disk Space Requirements

| Files | Grid 80³ | Grid 100³ | Grid 120³ |
|-------|----------|-----------|-----------|
| Per state (all 3 types) | ~150 MB | ~300 MB | ~500 MB |
| HOMO/LUMO (4 files) | ~200 MB | ~400 MB | ~700 MB |
| NTO molden (per state) | ~5 MB | ~5 MB | ~5 MB |

**Example:** 
- 3 states, all cube types, HOMO/LUMO, 80³ grid: ~650 MB
- 5 states, all cube types, HOMO/LUMO, 100³ grid: ~1.9 GB

---

## Troubleshooting

### Memory Issues

If you run out of memory:
1. Reduce `GRID_RESOLUTION` (e.g., from 100 to 80)
2. Reduce `NUM_EXCITED_STATES`
3. Reduce `STATES_TO_OUTPUT`
4. Disable some cube file types

### Slow Calculation

If calculation is too slow:
1. Increase `NUM_THREADS`
2. Reduce `NUM_EXCITED_STATES`
3. Use smaller basis set (6-31g instead of 6-31g*)
4. Reduce grid resolution

### File Not Found

If XYZ file not found:
1. Check `XYZ_FILE` path is correct
2. Use absolute path: `XYZ_FILE = '/full/path/to/PTCDA.xyz'`
3. Or ensure file is in same directory as script

---

## Workflow Recommendation

### Step 1: Test with H2O
```python
USE_XYZ = False
NUM_EXCITED_STATES = 5
STATES_TO_OUTPUT = [0]
GRID_RESOLUTION = [60, 60, 60]
```
Run time: ~2 minutes

### Step 2: Quick PTCDA Test
```python
USE_XYZ = True
NUM_EXCITED_STATES = 5
STATES_TO_OUTPUT = [0]
GRID_RESOLUTION = [60, 60, 60]
```
Run time: ~10 minutes

### Step 3: Full PTCDA Calculation
```python
USE_XYZ = True
NUM_EXCITED_STATES = 10
STATES_TO_OUTPUT = [0, 1, 2]
GRID_RESOLUTION = [80, 80, 80]
```
Run time: ~30 minutes

### Step 4: High Quality (if needed)
```python
BASIS_SET = '6-31g*'
GRID_RESOLUTION = [100, 100, 60]
```
Run time: ~60 minutes

---

## Quick Reference

**Most common settings for PTCDA:**

```python
ENABLE_PARALLEL = True
NUM_THREADS = 4
USE_XYZ = True
XYZ_FILE = 'PTCDA.xyz'
BASIS_SET = '6-31g'
NUM_EXCITED_STATES = 10
STATES_TO_OUTPUT = [0, 1, 2]
GENERATE_TRANSITION_DENSITY = True
GENERATE_EXCITED_DENSITY = True
GENERATE_DENSITY_DIFFERENCE = True
GENERATE_HOMO_LUMO = True
USE_GRID_RESOLUTION = True
GRID_RESOLUTION = [80, 80, 80]
ENABLE_NTO_ANALYSIS = True
NTO_STATES = [0, 1, 2]
OUTPUT_DIR = 'output'
```

Just modify the values at the top of the script and run!
