# Guide: Comparing HOMO×LUMO Multiplication vs Exact Methods

## Overview

This guide explains how to use `compare_homo_lumo_multiplication.py` to compare three methods of calculating transition density:

1. **Approximate (Multiplication):** ρ_approx(r) = 2 × φ_HOMO(r) × φ_LUMO(r)
2. **Exact (TDDFT):** ρ_exact(r) from TDDFT transition density matrix
3. **Exact (Analytical):** ρ_analytical(r) from HOMO⊗LUMO matrix evaluated properly

---

## Quick Start

### Step 1: Generate Cube Files

First, run the main script to generate the required cube files:

```bash
python tdm_calc_accurate.py
```

Make sure these settings are enabled:
```python
GENERATE_HOMO_LUMO = True
GENERATE_TRANSITION_DENSITY = True
STATES_TO_OUTPUT = [0]  # Include first state
```

This generates:
- `output/HOMO.cube`
- `output/LUMO.cube`
- `output/transition_density_state1.cube`
- `output/transition_HOMO_LUMO_analytical.cube`

### Step 2: Run Comparison Script

```bash
python compare_homo_lumo_multiplication.py
```

This will:
1. Load HOMO and LUMO cube files
2. Multiply them point-by-point: ρ(r) = 2 × φ_HOMO(r) × φ_LUMO(r)
3. Compare with TDDFT and analytical results
4. Generate statistics and difference cube files

### Step 3: Analyze Results

Check the output for similarity metrics:
```
Similarity (Multiplied vs TDDFT): 0.XXXX
```

- **> 0.95:** Approximation is excellent
- **0.85-0.95:** Approximation is good
- **0.70-0.85:** Approximation is fair
- **< 0.70:** Approximation is poor (basis overlap effects significant)

---

## What the Script Does

### 1. Reads Cube Files

Loads the volumetric data from:
- `HOMO.cube` → φ_HOMO(r) on 3D grid
- `LUMO.cube` → φ_LUMO(r) on 3D grid
- `transition_density_state1.cube` → ρ_TDDFT(r) from TDDFT
- `transition_HOMO_LUMO_analytical.cube` → ρ_analytical(r) from matrix

### 2. Performs Point-wise Multiplication

For each grid point r = (x, y, z):

```python
ρ_approx(r) = 2 × φ_HOMO(r) × φ_LUMO(r)
```

**Why the factor of 2?**

The transition density is symmetric:
```
T = |HOMO⟩⟨LUMO| + |LUMO⟩⟨HOMO|
```

In real space:
```
ρ_trans(r) = φ_HOMO(r) × φ_LUMO(r) + φ_LUMO(r) × φ_HOMO(r)
           = 2 × φ_HOMO(r) × φ_LUMO(r)
```

### 3. Calculates Statistics

For each comparison, computes:

**Basic statistics:**
- Min, max, mean, standard deviation
- Norm (Frobenius norm)

**Difference metrics:**
- Absolute difference: Δ = data1 - data2
- RMS difference: √(mean(Δ²))
- Relative difference: |Δ / data2|

**Similarity metrics:**
- Cosine similarity: ⟨data1, data2⟩ / (||data1|| × ||data2||)
- Correlation coefficient
- R² score (coefficient of determination)

### 4. Generates Output Files

**Main output:**
- `transition_HOMO_LUMO_multiplied.cube` - Result of multiplication

**Difference files (if enabled):**
- `difference_multiplied_minus_tddft.cube` - Shows approximation error
- `difference_multiplied_minus_analytical.cube` - Error vs analytical
- `difference_tddft_minus_analytical.cube` - Should be ~0 (both exact)

---

## Scientific Interpretation

### Why Compare These Methods?

**Question:** Is ρ_trans(r) ≈ φ_HOMO(r) × φ_LUMO(r) a good approximation?

**Answer:** It depends on the basis set and molecular system!

### The Approximation

**Approximate formula:**
```
ρ_approx(r) = 2 × φ_HOMO(r) × φ_LUMO(r)
```

**Exact formula:**
```
ρ_exact(r) = Σ_μν T_μν × χ_μ(r) × χ_ν(r)
```

where T_μν = C_μ^HOMO C_ν^LUMO + C_μ^LUMO C_ν^HOMO

### When Does the Approximation Work?

**Works well when:**
1. Basis functions are nearly orthogonal
2. Basis set is large and complete
3. HOMO and LUMO are well-separated spatially

**Works poorly when:**
1. Basis functions have significant overlap
2. Small basis set (e.g., STO-3G)
3. HOMO and LUMO are spatially close

### Expected Results

**For typical systems (6-31g or better basis):**
- Similarity: 0.85 - 0.95
- Correlation: > 0.90
- R²: > 0.80

**Interpretation:**
- The approximation captures the **qualitative** features correctly
- The **quantitative** values differ due to basis overlap
- The exact methods (TDDFT and Analytical) should agree to > 0.99

---

## Output Interpretation

### Example Output

```
======================================================================
STATISTICAL COMPARISON: Multiplied (approx) vs TDDFT (exact)
======================================================================

Multiplied (approx):
  Min:  -1.234567e-03
  Max:   2.345678e-03
  Mean:  1.234567e-06
  Std:   3.456789e-04
  Norm:  4.567890e-01

TDDFT (exact):
  Min:  -1.345678e-03
  Max:   2.456789e-03
  Mean:  1.345678e-06
  Std:   3.567890e-04
  Norm:  4.678901e-01

Difference (Multiplied - TDDFT):
  Min:  -2.345678e-04
  Max:   1.234567e-04
  Mean:  5.678901e-08
  Std:   4.567890e-05
  RMS:   4.567890e-05

Similarity Metrics:
  Cosine similarity: 0.923456
  Correlation coef:  0.945678
  R² score:          0.876543
```

### What to Look For

#### 1. Cosine Similarity

**Value:** 0.0 to 1.0

- **> 0.95:** Excellent agreement, approximation is very accurate
- **0.85-0.95:** Good agreement, approximation captures main features
- **0.70-0.85:** Fair agreement, some discrepancies
- **< 0.70:** Poor agreement, approximation not reliable

#### 2. Correlation Coefficient

**Value:** -1.0 to 1.0 (usually positive)

- **> 0.95:** Very strong linear relationship
- **0.85-0.95:** Strong linear relationship
- **< 0.85:** Weaker relationship

#### 3. R² Score

**Value:** 0.0 to 1.0

- **> 0.90:** Approximation explains >90% of variance
- **0.70-0.90:** Approximation explains 70-90% of variance
- **< 0.70:** Significant unexplained variance

#### 4. RMS Difference

**Value:** Depends on density magnitude

- Compare to the norm of the data
- RMS / Norm < 0.1 → Good agreement
- RMS / Norm > 0.2 → Significant differences

---

## Visualization

### Load in VMD

```bash
vmd output/HOMO.cube \
    output/LUMO.cube \
    output/transition_HOMO_LUMO_multiplied.cube \
    output/transition_density_state1.cube \
    output/transition_HOMO_LUMO_analytical.cube
```

### Create Isosurfaces

For each molecule:
```tcl
# Positive isosurface (red)
mol representation Isosurface 0.002 0 0 0 1 1
mol color ColorID 1
mol addrep 0

# Negative isosurface (blue)
mol representation Isosurface -0.002 0 0 0 1 1
mol color ColorID 0
mol addrep 0
```

### Compare Visually

1. **HOMO and LUMO:** See the frontier orbitals
2. **Multiplied:** Approximate transition density
3. **TDDFT:** Exact transition density from TDDFT
4. **Analytical:** Exact transition density from HOMO⊗LUMO matrix

**They should look similar but not identical!**

### View Differences

Load difference cube files:
```bash
vmd output/difference_multiplied_minus_tddft.cube
```

This shows where the approximation fails:
- Large positive/negative regions → Significant basis overlap effects
- Small uniform values → Approximation works well

---

## Configuration Options

Edit the script to customize:

```python
# Directory containing cube files
CUBE_DIR = 'output'

# Input files (automatically set based on CUBE_DIR)
HOMO_CUBE = os.path.join(CUBE_DIR, 'HOMO.cube')
LUMO_CUBE = os.path.join(CUBE_DIR, 'LUMO.cube')
TDDFT_TRANSITION_CUBE = os.path.join(CUBE_DIR, 'transition_density_state1.cube')
ANALYTICAL_CUBE = os.path.join(CUBE_DIR, 'transition_HOMO_LUMO_analytical.cube')

# Output file
OUTPUT_CUBE = os.path.join(CUBE_DIR, 'transition_HOMO_LUMO_multiplied.cube')

# Options
GENERATE_DIFFERENCE_CUBES = True  # Generate difference files
CALCULATE_STATISTICS = True       # Calculate detailed statistics
```

---

## Advanced Usage

### Compare Different States

To compare S₂ instead of S₁:

```python
TDDFT_TRANSITION_CUBE = os.path.join(CUBE_DIR, 'transition_density_state2.cube')
```

But note: S₂ might not be a simple HOMO→LUMO transition!

### Compare HOMO-1 → LUMO

```python
HOMO_CUBE = os.path.join(CUBE_DIR, 'HOMO-1.cube')
LUMO_CUBE = os.path.join(CUBE_DIR, 'LUMO.cube')
```

Then compare with the appropriate excited state.

### Use Different Basis Sets

Run `tdm_calc_accurate.py` with different `BASIS_SET` values:
- `'sto-3g'` - Minimal basis (poor approximation expected)
- `'6-31g'` - Small basis (fair approximation)
- `'6-31g*'` - Medium basis (good approximation)
- `'6-311g**'` - Large basis (excellent approximation)

Then compare the similarity metrics to see how basis set affects the approximation.

---

## Troubleshooting

### Error: "Required file not found"

**Solution:** Run `tdm_calc_accurate.py` first with:
```python
GENERATE_HOMO_LUMO = True
GENERATE_TRANSITION_DENSITY = True
STATES_TO_OUTPUT = [0]
```

### Error: "Grids don't match"

**Cause:** Cube files generated with different grid settings.

**Solution:** Regenerate all cube files with the same grid settings in `tdm_calc_accurate.py`.

### Low Similarity (<0.7)

**Possible causes:**
1. Small basis set (try larger basis)
2. Multi-configurational state (S₁ is not pure HOMO→LUMO)
3. Numerical issues

**Check:**
- Run NTO analysis to verify S₁ character
- Check HOMO→LUMO weight from verification output
- Try different basis set

### Very High Similarity (>0.99)

**Interpretation:** 
- The approximation works extremely well
- Basis functions are nearly orthogonal
- HOMO and LUMO are well-separated

This is actually good news!

---

## Scientific Notes

### Why the Approximation Differs

The exact transition density is:
```
ρ_exact(r) = Σ_μν T_μν χ_μ(r) χ_ν(r)
```

Expanding T_μν:
```
ρ_exact(r) = Σ_μν (C_μ^HOMO C_ν^LUMO + C_μ^LUMO C_ν^HOMO) χ_μ(r) χ_ν(r)
```

The approximation assumes:
```
Σ_μν C_μ^HOMO C_ν^LUMO χ_μ(r) χ_ν(r) ≈ [Σ_μ C_μ^HOMO χ_μ(r)] × [Σ_ν C_ν^LUMO χ_ν(r)]
                                        = φ_HOMO(r) × φ_LUMO(r)
```

This is only exact if χ_μ(r) χ_ν(r) = χ_μ(r) × χ_ν(r), which requires orthogonality.

### Basis Function Overlap

Gaussian basis functions are **not orthogonal**:
```
⟨χ_μ | χ_ν⟩ = S_μν ≠ δ_μν
```

The overlap matrix S_μν causes the approximation to differ from the exact result.

### Completeness

Even with orthogonal basis functions, a finite basis set is incomplete:
```
Σ_μ |χ_μ⟩⟨χ_μ| ≠ 1
```

This also contributes to the difference.

---

## Example Results

### Case 1: Good Approximation (6-31g basis)

```
Similarity: 0.923
Correlation: 0.945
R²: 0.876
RMS: 4.5e-05
```

**Interpretation:** Approximation captures 92% of the exact result. Good for qualitative analysis.

### Case 2: Excellent Approximation (6-311g** basis)

```
Similarity: 0.967
Correlation: 0.982
R²: 0.935
RMS: 2.1e-05
```

**Interpretation:** Approximation is very accurate. Larger basis set reduces overlap effects.

### Case 3: Poor Approximation (STO-3G basis)

```
Similarity: 0.734
Correlation: 0.812
R²: 0.659
RMS: 8.9e-05
```

**Interpretation:** Minimal basis set has significant overlap. Approximation not reliable.

---

## Summary

### Key Points

1. **The approximation ρ ≈ 2×φ_HOMO×φ_LUMO is qualitatively correct**
2. **Quantitative accuracy depends on basis set**
3. **The exact methods (TDDFT and Analytical) should agree perfectly**
4. **Use this script to understand basis set effects**

### Workflow

```
1. Run tdm_calc_accurate.py
   ↓
2. Run compare_homo_lumo_multiplication.py
   ↓
3. Check similarity metrics
   ↓
4. Visualize in VMD/Jmol
   ↓
5. Understand approximation quality
```

### Files Generated

- `transition_HOMO_LUMO_multiplied.cube` - Approximate method
- `difference_multiplied_minus_tddft.cube` - Error vs TDDFT
- `difference_multiplied_minus_analytical.cube` - Error vs analytical
- `difference_tddft_minus_analytical.cube` - Verification (should be ~0)

### Recommended Settings

For best results:
- Use 6-31g* or larger basis set
- Generate all cube files with same grid
- Enable difference cube file generation
- Visualize all files in VMD

---

## References

This script demonstrates concepts from:
- Quantum chemistry textbooks (Szabo & Ostlund, etc.)
- PySCF documentation
- Gaussian cube file format specification

The comparison helps understand the difference between approximate and exact quantum chemical calculations!
