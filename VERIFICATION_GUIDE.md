# Guide: Verifying Transition Density Matrix with HOMO/LUMO

## Overview

This guide explains how to analytically verify that the first excited state (S₁) transition density matrix corresponds to a HOMO→LUMO transition.

---

## Theory

### What is a HOMO→LUMO Transition?

For a simple single-electron excitation from HOMO to LUMO, the transition density matrix is:

**T = |HOMO⟩⟨LUMO| + |LUMO⟩⟨HOMO|**

In the atomic orbital (AO) basis:

```
T_μν = C_μ^HOMO × C_ν^LUMO + C_μ^LUMO × C_ν^HOMO
```

where:
- `C_μ^HOMO` = HOMO orbital coefficients in AO basis
- `C_ν^LUMO` = LUMO orbital coefficients in AO basis

### Why Verify?

1. **Sanity check:** Ensure TDDFT calculation is correct
2. **Understand excitation character:** Is S₁ really HOMO→LUMO?
3. **Multi-configurational states:** Some states mix multiple transitions
4. **Publication:** Reviewers often ask about excitation character

---

## Verification Methods

### Method 1: Visual Comparison (Qualitative)

**Files to compare:**
- `HOMO.cube` - HOMO orbital
- `LUMO.cube` - LUMO orbital
- `transition_density_state1.cube` - TDDFT transition density
- `transition_HOMO_LUMO_analytical.cube` - Analytical HOMO→LUMO (NEW!)

**What to look for:**

1. **HOMO orbital:**
   - Shows where electrons come FROM
   - These regions should correspond to NEGATIVE parts of transition density

2. **LUMO orbital:**
   - Shows where electrons go TO
   - These regions should correspond to POSITIVE parts of transition density

3. **Transition density:**
   - NEGATIVE regions = electron depletion (from HOMO)
   - POSITIVE regions = electron accumulation (to LUMO)

**In VMD:**
```tcl
# Load files
mol new output/HOMO.cube
mol new output/LUMO.cube
mol new output/transition_density_state1.cube
mol new output/transition_HOMO_LUMO_analytical.cube

# For each molecule, create isosurfaces
# Compare transition_density_state1 with transition_HOMO_LUMO_analytical
# They should look very similar if S1 is HOMO→LUMO
```

**In Jmol:**
```
load output/transition_density_state1.cube
isosurface ID "surf1" cutoff  0.002 color red
isosurface ID "surf2" cutoff -0.002 color blue

# Then load analytical version and compare
load output/transition_HOMO_LUMO_analytical.cube
isosurface ID "surf3" cutoff  0.002 color red
isosurface ID "surf4" cutoff -0.002 color blue
```

---

### Method 2: Quantitative Analysis (Analytical)

The script now automatically performs quantitative verification!

**Output example:**
```
----------------------------------------------------------------------
QUANTITATIVE VERIFICATION: S1 vs HOMO→LUMO
----------------------------------------------------------------------
Similarity (cosine): 0.9823
HOMO→LUMO weight: 0.8956 (89.6%)
HOMO→LUMO amplitude: 0.9463
✓ S1 is STRONGLY dominated by HOMO→LUMO transition
----------------------------------------------------------------------
```

**Metrics explained:**

#### 1. Similarity (Cosine Similarity)

**Formula:**
```
similarity = ⟨T_TDDFT | T_HOMO→LUMO⟩ / (||T_TDDFT|| × ||T_HOMO→LUMO||)
```

**Interpretation:**
- `1.0` = Perfect match (identical)
- `> 0.95` = Very strong HOMO→LUMO character
- `0.85-0.95` = Strong HOMO→LUMO character
- `0.70-0.85` = Significant HOMO→LUMO character, but mixed
- `< 0.70` = Multi-configurational (not pure HOMO→LUMO)

**What it measures:** How similar the TDDFT transition density is to the pure HOMO→LUMO transition density.

#### 2. HOMO→LUMO Weight

**Formula:**
```
weight = |X_HOMO→LUMO + Y_HOMO→LUMO|² / ||X + Y||²
```

where X and Y are the TDDFT amplitudes.

**Interpretation:**
- `> 0.8` = Strongly dominated by HOMO→LUMO
- `0.6-0.8` = Mostly HOMO→LUMO
- `0.4-0.6` = Significant but mixed
- `< 0.4` = Multi-configurational

**What it measures:** The fraction of the excitation that comes from HOMO→LUMO transition.

#### 3. HOMO→LUMO Amplitude

**Formula:**
```
amplitude = |X_HOMO→LUMO + Y_HOMO→LUMO|
```

**Interpretation:**
- `> 0.9` = Very strong single-excitation character
- `0.7-0.9` = Strong single-excitation
- `< 0.7` = Mixed character

**What it measures:** The raw amplitude of the HOMO→LUMO contribution.

---

## Interpretation Guide

### Case 1: Pure HOMO→LUMO Transition

**Metrics:**
```
Similarity: > 0.95
HOMO→LUMO weight: > 0.80
```

**Interpretation:**
- S₁ is essentially a pure HOMO→LUMO transition
- Single-configurational state
- Easy to interpret and visualize

**Example molecules:**
- Small conjugated systems (benzene, naphthalene)
- Simple π→π* transitions

### Case 2: Mostly HOMO→LUMO

**Metrics:**
```
Similarity: 0.85-0.95
HOMO→LUMO weight: 0.60-0.80
```

**Interpretation:**
- S₁ is dominated by HOMO→LUMO but has some mixing
- Minor contributions from other transitions (e.g., HOMO-1→LUMO, HOMO→LUMO+1)
- Still reasonable to describe as "HOMO→LUMO transition"

**Example molecules:**
- Medium-sized conjugated systems
- Some charge-transfer states

### Case 3: Mixed Character

**Metrics:**
```
Similarity: 0.70-0.85
HOMO→LUMO weight: 0.40-0.60
```

**Interpretation:**
- S₁ has significant HOMO→LUMO character but is mixed
- Multiple orbital pairs contribute
- Need to check NTO analysis for full picture

**Example molecules:**
- Large conjugated systems (like PTCDA)
- Degenerate or near-degenerate states

### Case 4: Multi-configurational

**Metrics:**
```
Similarity: < 0.70
HOMO→LUMO weight: < 0.40
```

**Interpretation:**
- S₁ is NOT a simple HOMO→LUMO transition
- Multiple configurations contribute significantly
- Cannot be described by single orbital pair
- **Must use NTO analysis**

**Example molecules:**
- Molecules with near-degenerate orbitals
- Charge-transfer states in large systems
- States with double-excitation character

---

## Step-by-Step Verification Procedure

### Step 1: Run the Script

```bash
python tdm_calc_accurate.py
```

Make sure:
- `GENERATE_HOMO_LUMO = True`
- `GENERATE_TRANSITION_DENSITY = True`
- `STATES_TO_OUTPUT` includes `[0]` (first state)

### Step 2: Check Quantitative Metrics

Look at the output:
```
QUANTITATIVE VERIFICATION: S1 vs HOMO→LUMO
Similarity (cosine): X.XXXX
HOMO→LUMO weight: X.XXXX (XX.X%)
```

**Decision tree:**
- Similarity > 0.95 AND weight > 0.80 → Pure HOMO→LUMO ✓
- Similarity > 0.85 AND weight > 0.60 → Mostly HOMO→LUMO ✓
- Similarity > 0.70 AND weight > 0.40 → Mixed character ⚠
- Otherwise → Multi-configurational ⚠ (check NTO)

### Step 3: Visual Verification

Load in VMD or Jmol:
1. `transition_density_state1.cube` (TDDFT result)
2. `transition_HOMO_LUMO_analytical.cube` (analytical HOMO→LUMO)

**They should look very similar if S₁ is HOMO→LUMO.**

### Step 4: Check NTO Analysis

If the state is mixed, check NTO weights:
```
State 1 (X.XX eV):
** Transition 1 **
  occ-alpha: 0.8956  # Should be > 0.8 for single-excitation
```

High NTO weight (> 0.8) confirms single-excitation character.

---

## Files Generated

### New Files for Verification:

1. **`HOMO.cube`** - HOMO orbital
2. **`LUMO.cube`** - LUMO orbital
3. **`HOMO-1.cube`** - Second highest occupied
4. **`LUMO+1.cube`** - Second lowest unoccupied
5. **`transition_HOMO_LUMO_analytical.cube`** ⭐ **NEW!**
   - Analytical HOMO→LUMO transition density
   - Calculated as: T = |HOMO⟩⟨LUMO| + |LUMO⟩⟨HOMO|
   - Compare directly with `transition_density_state1.cube`

---

## Mathematical Details

### Transition Density Matrix Construction

**From TDDFT (full TDDFT, not TDA):**
```python
X, Y = td.xy[state_id]  # TDDFT amplitudes
T_TDDFT = C_occ @ (X + Y) @ C_vir.T + C_vir @ (X + Y).T @ C_occ.T
```

**Analytical HOMO→LUMO:**
```python
T_analytical = np.outer(C_HOMO, C_LUMO) + np.outer(C_LUMO, C_HOMO)
```

where:
- `C_HOMO` = HOMO orbital coefficients (column vector)
- `C_LUMO` = LUMO orbital coefficients (column vector)
- `np.outer(a, b)` = outer product = a ⊗ b

### Similarity Calculation

**Frobenius inner product:**
```python
overlap = np.sum(T_TDDFT * T_analytical)  # Element-wise product, then sum
```

**Norms:**
```python
norm_TDDFT = np.linalg.norm(T_TDDFT)  # Frobenius norm
norm_analytical = np.linalg.norm(T_analytical)
```

**Cosine similarity:**
```python
similarity = overlap / (norm_TDDFT * norm_analytical)
```

This is equivalent to the cosine of the angle between the two matrices treated as vectors.

### Weight Calculation

**From TDDFT amplitudes:**
```python
X, Y = td.xy[0]  # First excited state
nocc, nvir = X.shape

# HOMO is last occupied (index nocc-1)
# LUMO is first virtual (index 0)
amplitude_HOMO_LUMO = abs(X[nocc-1, 0] + Y[nocc-1, 0])

# Total amplitude
total_amplitude = np.linalg.norm(X + Y)

# Weight (normalized squared amplitude)
weight = (amplitude_HOMO_LUMO / total_amplitude)**2
```

---

## Common Issues and Solutions

### Issue 1: Low Similarity but High Weight

**Symptoms:**
```
Similarity: 0.75
HOMO→LUMO weight: 0.85
```

**Cause:** The TDDFT transition density has the right orbital character but different spatial distribution due to:
- Orbital relaxation effects
- Correlation effects
- Basis set limitations

**Solution:** This is normal. The weight is more reliable for orbital character.

### Issue 2: High Similarity but Low Weight

**Symptoms:**
```
Similarity: 0.90
HOMO→LUMO weight: 0.55
```

**Cause:** Multiple transitions contribute with similar spatial patterns.

**Solution:** Check NTO analysis to see all contributing orbital pairs.

### Issue 3: Both Metrics Low

**Symptoms:**
```
Similarity: 0.60
HOMO→LUMO weight: 0.40
```

**Cause:** S₁ is genuinely multi-configurational.

**Solution:**
1. Check NTO analysis for dominant orbital pairs
2. Consider if this is a charge-transfer state
3. Try different XC functional (e.g., CAM-B3LYP)

---

## Example: PTCDA

For PTCDA, you might see:

**Scenario A: Localized π→π* transition**
```
Similarity: 0.92
HOMO→LUMO weight: 0.78
→ S1 is mostly HOMO→LUMO, good agreement
```

**Scenario B: Delocalized excitation**
```
Similarity: 0.75
HOMO→LUMO weight: 0.52
→ S1 has mixed character, check NTO
```

**What to do:**
1. Look at NTO weights
2. Visualize all cube files
3. Compare with experimental absorption spectrum
4. Check literature for PTCDA excited states

---

## Advanced: Manual Calculation in Python

If you want to verify manually:

```python
import numpy as np
from pyscf import gto, dft, tddft

# After running TDDFT...

# Get HOMO and LUMO
mo_occ = mf.mo_occ
homo_idx = np.where(mo_occ > 0)[0][-1]
lumo_idx = np.where(mo_occ == 0)[0][0]

# Get MO coefficients
C_homo = mf.mo_coeff[:, homo_idx]
C_lumo = mf.mo_coeff[:, lumo_idx]

# Construct analytical transition density
T_analytical = np.outer(C_homo, C_lumo) + np.outer(C_lumo, C_homo)

# Get TDDFT transition density for S1
X, Y = td.xy[0]
mo_coeff = mf.mo_coeff
nocc = X.shape[0]
nvir = X.shape[1]

t_dm1_mo = np.zeros((mo_coeff.shape[1], mo_coeff.shape[1]))
t_dm1_mo[:nocc, nocc:] = (X + Y).reshape(nocc, nvir)
t_dm1_mo[nocc:, :nocc] = (X + Y).reshape(nocc, nvir).T
T_TDDFT = mo_coeff @ t_dm1_mo @ mo_coeff.T

# Calculate similarity
overlap = np.sum(T_TDDFT * T_analytical)
similarity = overlap / (np.linalg.norm(T_TDDFT) * np.linalg.norm(T_analytical))

print(f"Similarity: {similarity:.4f}")
```

---

## Summary

### Quick Checklist:

- [ ] Run script with `GENERATE_HOMO_LUMO = True`
- [ ] Check quantitative metrics in output
- [ ] Compare `transition_density_state1.cube` with `transition_HOMO_LUMO_analytical.cube`
- [ ] If similarity > 0.95 and weight > 0.80: ✓ Pure HOMO→LUMO
- [ ] If metrics lower: Check NTO analysis
- [ ] Visualize all files in VMD/Jmol

### Key Files:

1. `transition_density_state1.cube` - TDDFT result
2. `transition_HOMO_LUMO_analytical.cube` - Analytical HOMO→LUMO
3. `HOMO.cube`, `LUMO.cube` - Frontier orbitals
4. `nto_state_1.molden` - NTO analysis

### Interpretation:

| Similarity | Weight | Interpretation |
|------------|--------|----------------|
| > 0.95 | > 0.80 | Pure HOMO→LUMO ✓ |
| 0.85-0.95 | 0.60-0.80 | Mostly HOMO→LUMO ✓ |
| 0.70-0.85 | 0.40-0.60 | Mixed character ⚠ |
| < 0.70 | < 0.40 | Multi-configurational ⚠ |

The script now does all of this automatically! Just run it and check the output.
