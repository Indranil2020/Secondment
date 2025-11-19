# Transition Contribution Analysis Tool

## Overview

Analyze which orbital pairs contribute to each excited state in TDDFT calculations - similar to GPAW's transition analysis!

**Key Question:** Is S₁ really 100% HOMO→LUMO, or are there other contributions?

---

## Quick Start

```bash
python analyze_transition_contributions.py
```

**Output:**
```
STATE 1: 2.8456 eV
Rank   Transition           Weight       Percentage   Cumulative  
----------------------------------------------------------------------
1      HOMO → LUMO          0.856234     85.62%       85.62%
2      HOMO-1 → LUMO        0.089123     8.91%        94.53%
3      HOMO → LUMO+1        0.034567     3.46%        97.99%
```

---

## What It Does

### 1. Analyzes Contributions

For each excited state, shows:
- Which orbital pairs (i→a) contribute
- Percentage contribution of each pair
- Cumulative contributions

### 2. Generates Cube Files

Creates transition density maps for individual orbital pairs:
- `transition_pair_state1_HOMO_to_LUMO.cube`
- `transition_pair_state1_HOMOm1_to_LUMO.cube`
- etc.

### 3. Saves Tables

Exports detailed contribution tables to `contribution_tables.txt`

---

## Configuration

Edit the script:

```python
# Which states to analyze
STATES_TO_ANALYZE = [0, 1, 2]  # First 3 states

# Show contributions > 1%
CONTRIBUTION_THRESHOLD = 0.01

# Generate cube files for top 3 pairs per state
MAX_PAIRS_PER_STATE = 3
PAIR_CONTRIBUTION_CUTOFF = 0.05  # Only if >5%
```

---

## Understanding Results

### Single-Configurational State ✓

```
1      HOMO → LUMO          0.856234     85.62%
2      HOMO-1 → LUMO        0.089123     8.91%
```

**Interpretation:**
- One dominant transition (>80%)
- Easy to interpret as "HOMO→LUMO transition"
- Other pairs contribute minimally

### Multi-Configurational State ⚠

```
1      HOMO → LUMO          0.456234     45.62%
2      HOMO-1 → LUMO        0.289123     28.91%
3      HOMO → LUMO+1        0.134567     13.46%
```

**Interpretation:**
- Multiple important transitions
- Cannot be described by single orbital pair
- Need to consider several contributions

---

## Visualization

### Load in VMD

```bash
vmd transition_analysis/transition_pair_state1_HOMO_to_LUMO.cube \
    transition_analysis/transition_pair_state1_HOMOm1_to_LUMO.cube \
    output/transition_density_state1.cube
```

### What to Compare

1. **Individual pairs** - See each orbital pair's spatial contribution
2. **Full TDDFT** - See total transition density
3. **Verify** - Dominant pair should resemble full density

---

## Output Files

### 1. Contribution Tables

**File:** `transition_analysis/contribution_tables.txt`

**Content:**
```
STATE 1: 2.8456 eV
Rank   Transition           Weight       Percentage   Cumulative  
1      HOMO → LUMO          0.856234     85.62%       85.62%
...
```

### 2. Orbital Pair Cube Files

**Naming:** `transition_pair_state{N}_{OCC}_to_{VIR}.cube`

**Examples:**
- `transition_pair_state1_HOMO_to_LUMO.cube`
- `transition_pair_state1_HOMOm1_to_LUMO.cube`
- `transition_pair_state2_HOMO_to_LUMOp1.cube`

---

## Use Cases

### 1. Verify Transition Character

**Question:** Is S₁ pure HOMO→LUMO?

**Check:**
- If HOMO→LUMO > 80% → Yes, pure ✓
- If HOMO→LUMO < 60% → No, mixed ⚠

### 2. Identify Mixed States

**Question:** Which states have multi-configurational character?

**Look for:**
- No single dominant pair (<50%)
- Multiple significant pairs (>10% each)

### 3. Visualize Contributions

**Question:** Where does each orbital pair contribute spatially?

**Do:**
- Generate cube files for top pairs
- Load in VMD
- Compare spatial distributions

### 4. Publication

**Use contribution tables to report:**
```
S₁ (2.85 eV): HOMO→LUMO (85.6%), HOMO-1→LUMO (8.9%)
S₂ (3.12 eV): HOMO→LUMO+1 (45.2%), HOMO-1→LUMO (32.1%)
```

---

## How It Works

### TDDFT Amplitudes

TDDFT gives amplitudes X and Y for each orbital pair transition:

```
|Ψ_excited⟩ = Σ_(i,a) (X_(i,a) + Y_(i,a)) |i→a⟩
```

### Contribution Weights

The contribution of each pair is:

```
w_(i,a) = (X_(i,a) + Y_(i,a))²
```

Normalized so Σ w_(i,a) = 1

### Transition Density for Each Pair

For each orbital pair (i→a):

```
T^(i→a)_μν = C^i_μ × C^a_ν + C^a_μ × C^i_ν
```

### Full Transition Density

The full TDDFT transition density is:

```
T_total = Σ_(i,a) (X_(i,a) + Y_(i,a)) × T^(i→a)
```

---

## Comparison with GPAW

### GPAW Approach

```python
# In GPAW
calc.write('transitions.dat')
```

Shows orbital pair contributions in text format.

### This Script (PySCF)

**Equivalent functionality:**
- Analyzes TDDFT amplitudes
- Shows orbital pair contributions
- **Plus:** Generates cube files for visualization

**Advantages:**
- More flexible (any basis set)
- Visual output (cube files)
- Detailed tables

---

## Example: PTCDA

### Expected Results

**S₁:**
```
1      HOMO → LUMO          ~75-85%
2      HOMO-1 → LUMO        ~8-12%
3      HOMO → LUMO+1        ~3-6%
```

**Interpretation:** S₁ is mostly HOMO→LUMO (π→π*)

**S₂:**
```
1      HOMO-1 → LUMO        ~45-55%
2      HOMO → LUMO+1        ~25-35%
3      HOMO → LUMO          ~8-12%
```

**Interpretation:** S₂ has mixed character

---

## Tips

### For Quick Analysis

```python
STATES_TO_ANALYZE = [0]  # Just S₁
CONTRIBUTION_THRESHOLD = 0.05  # Only >5%
GENERATE_PAIR_CUBES = False  # No cubes
```

**Fast:** Just prints contribution table

### For Detailed Analysis

```python
STATES_TO_ANALYZE = range(5)  # First 5 states
CONTRIBUTION_THRESHOLD = 0.01  # Show >1%
GENERATE_PAIR_CUBES = True  # Generate cubes
MAX_PAIRS_PER_STATE = 5  # Top 5 pairs
```

**Thorough:** Full analysis with visualization

### For Publication

```python
STATES_TO_ANALYZE = range(10)  # All states
SAVE_CONTRIBUTION_TABLE = True  # Save tables
GENERATE_PAIR_CUBES = True  # For figures
MAX_PAIRS_PER_STATE = 3  # Top 3 for clarity
```

**Complete:** Everything for paper

---

## Troubleshooting

### "Contributions don't sum to 100%"

**Normal!** The table shows top N contributions. Check "Total weight analyzed" - should be ~1.0.

### "Too many cube files"

**Solution:**
```python
MAX_PAIRS_PER_STATE = 2  # Reduce
PAIR_CONTRIBUTION_CUTOFF = 0.10  # Increase to 10%
```

### "No significant contributions"

**Solution:**
```python
CONTRIBUTION_THRESHOLD = 0.001  # Lower threshold
```

---

## Summary

### What You Get

| Feature | Description |
|---------|-------------|
| **Contribution tables** | Which orbital pairs contribute |
| **Percentage weights** | How much each pair contributes |
| **Cube files** | Spatial distribution of each pair |
| **Verification** | Confirm transition character |

### Workflow

```
Configure → Run → Check tables → Visualize → Understand
```

### Key Insight

**Not all transitions are pure HOMO→LUMO!**

This tool shows you:
- ✅ Which orbital pairs contribute
- ✅ How much each contributes
- ✅ Where each contributes spatially

Perfect for understanding excited state character! 🎓

---

## Documentation

- **Quick start:** This file
- **Detailed guide:** `TRANSITION_ANALYSIS_GUIDE.md`
- **Theory:** See GPAW tutorial and PySCF examples

---

## Example Output

```
======================================================================
TRANSITION CONTRIBUTION ANALYSIS
======================================================================

STATE 1: 2.8456 eV
======================================================================
Rank   Transition           Weight       Percentage   Cumulative  
----------------------------------------------------------------------
1      HOMO → LUMO          0.856234     85.62%       85.62%
2      HOMO-1 → LUMO        0.089123     8.91%        94.53%
3      HOMO → LUMO+1        0.034567     3.46%        97.99%
----------------------------------------------------------------------

Generated files:
  ✓ transition_pair_state1_HOMO_to_LUMO.cube (85.6%)
  ✓ transition_pair_state1_HOMOm1_to_LUMO.cube (8.9%)

Key finding: S1 is STRONGLY dominated by HOMO→LUMO transition
```

This is exactly what you asked for - orbital contribution analysis like GPAW! 🚀
