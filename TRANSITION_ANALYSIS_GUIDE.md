# Guide: Transition Contribution Analysis

## Overview

This guide explains how to use `analyze_transition_contributions.py` to analyze which orbital pairs contribute to each excited state, similar to GPAW's transition contribution analysis.

---

## What This Script Does

### 1. **Analyzes Orbital Pair Contributions**

For each excited state, shows which orbital transitions (i→a) contribute:

```
STATE 1: 2.8456 eV
Rank   Transition           Weight       Percentage   Cumulative  
----------------------------------------------------------------------
1      HOMO → LUMO          0.856234     85.62%       85.62%
2      HOMO-1 → LUMO        0.089123     8.91%        94.53%
3      HOMO → LUMO+1        0.034567     3.46%        97.99%
...
```

### 2. **Generates Transition Density Maps for Individual Pairs**

Creates cube files showing the transition density for each dominant orbital pair:

- `transition_pair_state1_HOMO_to_LUMO.cube` - HOMO→LUMO contribution
- `transition_pair_state1_HOMOm1_to_LUMO.cube` - HOMO-1→LUMO contribution
- etc.

### 3. **Saves Detailed Tables**

Exports contribution tables to `contribution_tables.txt` for reference.

---

## Quick Start

### Step 1: Configure the Script

Edit `analyze_transition_contributions.py`:

```python
# Which states to analyze
STATES_TO_ANALYZE = [0, 1, 2]  # First 3 states (0-indexed)

# Analysis settings
CONTRIBUTION_THRESHOLD = 0.01  # Show contributions > 1%
TOP_N_CONTRIBUTIONS = 10       # Show top 10 orbital pairs

# Cube file generation
GENERATE_PAIR_CUBES = True     # Generate cube files
MAX_PAIRS_PER_STATE = 3        # Top 3 pairs per state
PAIR_CONTRIBUTION_CUTOFF = 0.05  # Only if contribution > 5%
```

### Step 2: Run the Script

```bash
python analyze_transition_contributions.py
```

### Step 3: Check the Output

**Console output:**
```
STATE 1: 2.8456 eV
Rank   Transition           Weight       Percentage   Cumulative  
1      HOMO → LUMO          0.856234     85.62%       85.62%
2      HOMO-1 → LUMO        0.089123     8.91%        94.53%
...
```

**Files generated:**
- `transition_analysis/contribution_tables.txt` - Detailed tables
- `transition_analysis/transition_pair_state1_HOMO_to_LUMO.cube` - Individual pair densities

---

## Understanding the Output

### Contribution Table

```
Rank   Transition           Weight       Percentage   Cumulative  
----------------------------------------------------------------------
1      HOMO → LUMO          0.856234     85.62%       85.62%
```

**Columns explained:**

- **Rank:** Importance ranking (1 = most important)
- **Transition:** Which orbital pair (e.g., HOMO → LUMO)
- **Weight:** Normalized squared amplitude (sums to 1.0)
- **Percentage:** Contribution percentage
- **Cumulative:** Running total

### Interpretation

**Single-configurational state:**
```
1      HOMO → LUMO          0.856234     85.62%       85.62%
2      HOMO-1 → LUMO        0.089123     8.91%        94.53%
```
- One dominant transition (>80%)
- Easy to interpret
- S₁ is "pure" HOMO→LUMO

**Multi-configurational state:**
```
1      HOMO → LUMO          0.456234     45.62%       45.62%
2      HOMO-1 → LUMO        0.289123     28.91%       74.53%
3      HOMO → LUMO+1        0.134567     13.46%       87.99%
```
- Multiple important transitions
- More complex character
- Need to consider several pairs

---

## Orbital Pair Transition Densities

### What Are They?

For each orbital pair (i→a), the script calculates:

```
T_μν^(i→a) = C_μ^i × C_ν^a + C_μ^a × C_ν^i
```

This shows the transition density **only for that specific pair**.

### How They Relate to Full Transition Density

The full TDDFT transition density is a **weighted sum** of all pair contributions:

```
T_total = Σ_(i,a) w_(i,a) × T^(i→a)
```

where w_(i,a) are the TDDFT amplitudes.

### Example

For S₁ with:
- HOMO→LUMO: 85.6%
- HOMO-1→LUMO: 8.9%
- HOMO→LUMO+1: 3.5%

The full transition density is approximately:

```
T_S1 ≈ 0.856 × T^(HOMO→LUMO) + 0.089 × T^(HOMO-1→LUMO) + 0.035 × T^(HOMO→LUMO+1)
```

---

## Visualization

### Load in VMD

```bash
# Load all contributions for State 1
vmd transition_analysis/transition_pair_state1_HOMO_to_LUMO.cube \
    transition_analysis/transition_pair_state1_HOMOm1_to_LUMO.cube \
    output/transition_density_state1.cube
```

### Compare

1. **Individual pairs:** See each orbital pair's contribution
2. **Full TDDFT:** See the total transition density
3. **Verify:** The dominant pair should look similar to the full density

### Create Isosurfaces

```tcl
# For each molecule
mol representation Isosurface 0.002 0 0 0 1 1
mol color ColorID 1
mol addrep 0

mol representation Isosurface -0.002 0 0 0 1 1
mol color ColorID 0
mol addrep 0
```

---

## Configuration Options

### Analysis Settings

```python
# Which states to analyze
STATES_TO_ANALYZE = [0, 1, 2]  # 0-indexed (State 1, 2, 3)

# Contribution threshold
CONTRIBUTION_THRESHOLD = 0.01  # Show pairs with >1% contribution

# Number of pairs to show
TOP_N_CONTRIBUTIONS = 10  # Show top 10 pairs
```

**Recommendations:**
- `CONTRIBUTION_THRESHOLD = 0.01` (1%) - Good for detailed analysis
- `CONTRIBUTION_THRESHOLD = 0.05` (5%) - Focus on important pairs
- `TOP_N_CONTRIBUTIONS = 10` - Usually sufficient

### Cube File Generation

```python
# Enable/disable cube file generation
GENERATE_PAIR_CUBES = True

# How many pairs per state
MAX_PAIRS_PER_STATE = 3  # Top 3 pairs

# Minimum contribution to generate cube
PAIR_CONTRIBUTION_CUTOFF = 0.05  # Only if >5%
```

**Recommendations:**
- `MAX_PAIRS_PER_STATE = 3` - Good balance
- `PAIR_CONTRIBUTION_CUTOFF = 0.05` - Avoid tiny contributions
- Set to `0.10` (10%) if you only want dominant pairs

### Grid Settings

```python
# Grid resolution
USE_GRID_RESOLUTION = False
GRID_RESOLUTION = [80, 80, 80]

# Or use box dimensions
BOX_MARGIN = 4.0  # Angstrom
GRID_SPACING = 0.2  # Angstrom
```

**Same as `tdm_calc_accurate.py`** - Use consistent settings!

---

## Use Cases

### Case 1: Verify HOMO→LUMO Character

**Question:** Is S₁ really a pure HOMO→LUMO transition?

**How to check:**
```python
STATES_TO_ANALYZE = [0]  # Just S₁
CONTRIBUTION_THRESHOLD = 0.01
```

**Look for:**
```
1      HOMO → LUMO          0.856234     85.62%
```

- **>80%:** Pure HOMO→LUMO ✓
- **60-80%:** Mostly HOMO→LUMO
- **<60%:** Multi-configurational

### Case 2: Identify Multi-configurational States

**Question:** Which states have mixed character?

**How to check:**
```python
STATES_TO_ANALYZE = range(10)  # All states
TOP_N_CONTRIBUTIONS = 5
```

**Look for:**
- States with no single dominant pair (<50%)
- States with multiple significant pairs (>10% each)

### Case 3: Visualize Individual Contributions

**Question:** What does each orbital pair contribute spatially?

**How to do it:**
```python
GENERATE_PAIR_CUBES = True
MAX_PAIRS_PER_STATE = 5
PAIR_CONTRIBUTION_CUTOFF = 0.05
```

**Then visualize in VMD:**
- See spatial distribution of each pair
- Compare with full transition density
- Understand where each pair contributes

### Case 4: Compare with Experimental Data

**Question:** Which orbital transitions explain the observed spectrum?

**How to use:**
1. Run analysis for all states
2. Check contribution tables
3. Match dominant transitions with experimental peaks
4. Assign spectral features to specific orbital pairs

---

## Output Files

### 1. `contribution_tables.txt`

**Content:**
```
======================================================================
ORBITAL PAIR CONTRIBUTIONS TO EXCITED STATES
======================================================================

======================================================================
STATE 1: 2.8456 eV
======================================================================
Rank   Transition           Weight       Percentage   Cumulative  
----------------------------------------------------------------------
1      HOMO → LUMO          0.856234     85.62%       85.62%
2      HOMO-1 → LUMO        0.089123     8.91%        94.53%
...
```

**Use for:**
- Reference
- Publications
- Further analysis

### 2. Orbital Pair Cube Files

**Naming convention:**
```
transition_pair_state{N}_{OCC}_to_{VIR}.cube
```

**Examples:**
- `transition_pair_state1_HOMO_to_LUMO.cube`
- `transition_pair_state1_HOMOm1_to_LUMO.cube` (HOMO-1)
- `transition_pair_state2_HOMO_to_LUMOp1.cube` (LUMO+1)

**Use for:**
- Visualization
- Spatial analysis
- Comparison with full TDDFT

---

## Advanced Usage

### Analyze Specific Orbital Pairs

To see contribution of a specific pair (e.g., HOMO-2 → LUMO+1):

1. Run the script
2. Check `contribution_tables.txt`
3. Look for the pair in the table
4. If significant (>5%), cube file will be generated

### Compare Different XC Functionals

Run with different functionals to see how contributions change:

```python
# Run 1
XC_FUNCTIONAL = 'b3lyp'

# Run 2
XC_FUNCTIONAL = 'cam-b3lyp'
```

Compare contribution tables to see functional dependence.

### Export for Publication

The contribution tables are formatted for easy inclusion in papers:

```
State 1 (2.85 eV): HOMO→LUMO (85.6%), HOMO-1→LUMO (8.9%)
State 2 (3.12 eV): HOMO→LUMO+1 (45.2%), HOMO-1→LUMO (32.1%)
```

---

## Comparison with GPAW

### GPAW Approach

GPAW's `lcaotddft` module provides:
```python
calc.write('transitions.dat')
```

Which shows orbital pair contributions.

### PySCF Approach (This Script)

**Equivalent functionality:**
- Analyzes TDDFT amplitudes (X, Y)
- Calculates squared amplitudes as weights
- Generates transition density for each pair

**Advantages:**
- More flexible (can analyze any state)
- Generates cube files for visualization
- Works with any basis set

**Differences:**
- GPAW uses LCAO basis (atomic orbitals)
- PySCF uses Gaussian basis (more flexible)
- Results should be very similar

---

## Troubleshooting

### Issue: "No significant contributions found"

**Cause:** Threshold too high

**Solution:**
```python
CONTRIBUTION_THRESHOLD = 0.001  # Lower to 0.1%
```

### Issue: "Too many cube files generated"

**Cause:** Too many pairs selected

**Solution:**
```python
MAX_PAIRS_PER_STATE = 2  # Reduce to top 2
PAIR_CONTRIBUTION_CUTOFF = 0.10  # Increase to 10%
```

### Issue: "Contributions don't sum to 100%"

**Cause:** Only showing top N contributions

**Explanation:** The table shows `TOP_N_CONTRIBUTIONS` pairs, but there may be more small contributions not shown. The "Total weight analyzed" shows the sum of all contributions (should be ~1.0).

### Issue: "Cube files look wrong"

**Cause:** Grid settings mismatch

**Solution:** Use same grid settings as `tdm_calc_accurate.py`:
```python
BOX_MARGIN = 4.0
GRID_SPACING = 0.2
```

---

## Example: PTCDA Analysis

### Expected Results

For PTCDA with B3LYP/6-31g:

**S₁ (first excited state):**
```
1      HOMO → LUMO          ~70-85%
2      HOMO-1 → LUMO        ~5-15%
3      HOMO → LUMO+1        ~2-8%
```

**S₂ (second excited state):**
```
1      HOMO-1 → LUMO        ~40-60%
2      HOMO → LUMO+1        ~20-40%
3      HOMO → LUMO          ~5-15%
```

**Interpretation:**
- S₁ is mostly HOMO→LUMO (π→π*)
- S₂ has mixed character (multi-configurational)

### Visualization

```bash
# Load S1 contributions
vmd transition_analysis/transition_pair_state1_HOMO_to_LUMO.cube \
    transition_analysis/transition_pair_state1_HOMOm1_to_LUMO.cube \
    output/transition_density_state1.cube
```

**What you'll see:**
- HOMO→LUMO dominates (looks similar to full density)
- HOMO-1→LUMO adds small corrections
- Sum approximately equals full TDDFT density

---

## Summary

### Key Features

1. ✅ **Analyzes orbital pair contributions** - Shows which transitions matter
2. ✅ **Generates individual pair densities** - Visualize each contribution
3. ✅ **Saves detailed tables** - For reference and publication
4. ✅ **Flexible configuration** - Customize analysis depth

### Workflow

```
1. Configure analysis settings
   ↓
2. Run script
   ↓
3. Check contribution tables
   ↓
4. Visualize dominant pairs
   ↓
5. Understand transition character
```

### Files Generated

- `contribution_tables.txt` - Detailed contribution data
- `transition_pair_state{N}_{OCC}_to_{VIR}.cube` - Individual pair densities

### Comparison with GPAW

**Similar to GPAW's transition analysis but:**
- More flexible (any basis set)
- Generates cube files
- Works with PySCF ecosystem

---

## References

- GPAW TDDFT tutorial: https://gpaw.readthedocs.io/tutorialsexercises/opticalresponse/tddft/lcaotddft.html
- PySCF TDDFT examples: https://github.com/pyscf/pyscf/tree/master/examples/tddft
- NTO analysis: Similar concept, but focuses on natural transition orbitals

This script provides the orbital pair analysis that GPAW users are familiar with, adapted for PySCF! 🎓
