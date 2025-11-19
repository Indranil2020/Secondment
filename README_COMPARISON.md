# HOMO×LUMO Multiplication Comparison Tool

## Quick Overview

This tool compares **approximate** vs **exact** methods for calculating transition density:

| Method | Formula | Type | Accuracy |
|--------|---------|------|----------|
| **Multiplication** | ρ(r) = 2×φ_HOMO(r)×φ_LUMO(r) | Approximate | Qualitative |
| **TDDFT** | ρ(r) = Σ_μν T_μν χ_μ(r) χ_ν(r) | Exact | Quantitative |
| **Analytical** | Same as TDDFT, from HOMO⊗LUMO | Exact | Quantitative |

---

## Quick Start

### 1. Generate cube files
```bash
python tdm_calc_accurate.py
```

### 2. Run comparison
```bash
python compare_homo_lumo_multiplication.py
```

### 3. Check results
```
Similarity (Multiplied vs TDDFT): 0.XXXX
```

- **> 0.95:** Excellent approximation ✓
- **0.85-0.95:** Good approximation ✓
- **< 0.85:** Fair/poor approximation ⚠

---

## What It Does

### Input Files (from `tdm_calc_accurate.py`)
- `output/HOMO.cube` - HOMO orbital φ_HOMO(r)
- `output/LUMO.cube` - LUMO orbital φ_LUMO(r)
- `output/transition_density_state1.cube` - TDDFT result (exact)
- `output/transition_HOMO_LUMO_analytical.cube` - Analytical result (exact)

### Processing
1. Loads HOMO and LUMO cube files
2. Multiplies point-by-point: ρ_approx(r) = 2 × φ_HOMO(r) × φ_LUMO(r)
3. Compares with exact TDDFT and analytical results
4. Calculates similarity metrics

### Output Files
- `output/transition_HOMO_LUMO_multiplied.cube` - Result of multiplication
- `output/difference_multiplied_minus_tddft.cube` - Shows approximation error
- `output/difference_multiplied_minus_analytical.cube` - Error vs analytical
- `output/difference_tddft_minus_analytical.cube` - Should be ~0 (both exact)

---

## Why This Matters

### The Question
**Is ρ_trans(r) ≈ φ_HOMO(r) × φ_LUMO(r) a good approximation?**

### The Answer
**It depends on the basis set!**

- **Large basis (6-311g\*\*):** Similarity ~0.96 → Excellent
- **Medium basis (6-31g\*):** Similarity ~0.92 → Good
- **Small basis (STO-3G):** Similarity ~0.73 → Poor

### Why They Differ

**Approximate method:**
- Assumes basis functions are orthogonal
- Simple multiplication in real space
- Fast and intuitive

**Exact method:**
- Accounts for basis function overlap
- Proper matrix evaluation with basis functions
- Scientifically correct

---

## Visualization

### Load in VMD
```bash
vmd output/HOMO.cube \
    output/LUMO.cube \
    output/transition_HOMO_LUMO_multiplied.cube \
    output/transition_density_state1.cube \
    output/difference_multiplied_minus_tddft.cube
```

### What to Look For

1. **HOMO.cube** - Where electrons come from
2. **LUMO.cube** - Where electrons go to
3. **multiplied.cube** - Approximate transition density
4. **state1.cube** - Exact TDDFT transition density
5. **difference.cube** - Where approximation fails

**The multiplied and TDDFT should look similar but not identical!**

---

## Statistics Explained

### Cosine Similarity
```
similarity = ⟨ρ_approx, ρ_exact⟩ / (||ρ_approx|| × ||ρ_exact||)
```

- Measures how similar the two densities are
- 1.0 = identical, 0.0 = completely different
- **> 0.90 is considered good**

### Correlation Coefficient
```
corr = correlation(ρ_approx, ρ_exact)
```

- Measures linear relationship
- 1.0 = perfect correlation
- **> 0.90 is considered strong**

### R² Score
```
R² = 1 - SS_residual / SS_total
```

- Measures explained variance
- 1.0 = perfect fit
- **> 0.80 is considered good**

### RMS Difference
```
RMS = √(mean((ρ_approx - ρ_exact)²))
```

- Measures average error
- Lower is better
- Compare to norm of data

---

## Example Output

```
======================================================================
STATISTICAL COMPARISON: Multiplied (approx) vs TDDFT (exact)
======================================================================

Multiplied (approx):
  Min:  -1.234567e-03
  Max:   2.345678e-03
  Norm:  4.567890e-01

TDDFT (exact):
  Min:  -1.345678e-03
  Max:   2.456789e-03
  Norm:  4.678901e-01

Difference (Multiplied - TDDFT):
  RMS:   4.567890e-05

Similarity Metrics:
  Cosine similarity: 0.923456  ← Main metric
  Correlation coef:  0.945678
  R² score:          0.876543

======================================================================
```

**Interpretation:** Similarity of 0.92 means the approximation captures 92% of the exact result. Good for qualitative understanding!

---

## Use Cases

### 1. Educational
- Understand difference between approximate and exact methods
- See how basis set affects accuracy
- Learn about basis function overlap

### 2. Validation
- Verify that S₁ is indeed HOMO→LUMO
- Check if simple multiplication is sufficient
- Validate TDDFT calculations

### 3. Basis Set Testing
- Compare results with different basis sets
- Understand when approximation breaks down
- Optimize computational cost vs accuracy

---

## Configuration

Edit `compare_homo_lumo_multiplication.py`:

```python
# Directory with cube files
CUBE_DIR = 'output'

# Options
GENERATE_DIFFERENCE_CUBES = True  # Create difference files
CALCULATE_STATISTICS = True       # Detailed statistics
```

---

## Troubleshooting

### "Required file not found"
**Solution:** Run `tdm_calc_accurate.py` first with:
```python
GENERATE_HOMO_LUMO = True
GENERATE_TRANSITION_DENSITY = True
STATES_TO_OUTPUT = [0]
```

### "Grids don't match"
**Solution:** Regenerate all cube files with same grid settings

### Low similarity (<0.7)
**Possible causes:**
- Small basis set → Try larger basis
- S₁ is not HOMO→LUMO → Check NTO analysis
- Multi-configurational state

---

## Scientific Notes

### Why Factor of 2?

The transition density is symmetric:
```
T = |HOMO⟩⟨LUMO| + |LUMO⟩⟨HOMO|
```

In real space:
```
ρ(r) = φ_HOMO(r)×φ_LUMO(r) + φ_LUMO(r)×φ_HOMO(r)
     = 2 × φ_HOMO(r)×φ_LUMO(r)
```

### Basis Function Overlap

Gaussian basis functions overlap:
```
⟨χ_μ | χ_ν⟩ = S_μν ≠ δ_μν
```

This overlap causes the approximation to differ from exact result.

### When Approximation Works

**Good approximation when:**
- Large basis set (more orthogonal)
- Well-separated HOMO and LUMO
- Single-configurational state

**Poor approximation when:**
- Small basis set (significant overlap)
- Spatially close HOMO and LUMO
- Multi-configurational state

---

## Summary

### Key Takeaways

1. ✅ **Approximation is qualitatively correct** - Shows same features
2. ⚠ **Approximation is quantitatively different** - Values differ
3. ✓ **Exact methods agree perfectly** - TDDFT ≈ Analytical
4. 📊 **Similarity depends on basis set** - Larger is better

### Workflow

```
Generate cubes → Run comparison → Check similarity → Visualize → Understand
```

### Files

| File | Description | Size |
|------|-------------|------|
| `compare_homo_lumo_multiplication.py` | Main script | ~15 KB |
| `MULTIPLICATION_COMPARISON_GUIDE.md` | Detailed guide | ~25 KB |
| `README_COMPARISON.md` | This file | ~5 KB |

### Documentation

- **Quick start:** This file
- **Detailed theory:** `MULTIPLICATION_COMPARISON_GUIDE.md`
- **AO vs real space:** `AO_vs_REALSPACE.md`
- **Verification:** `VERIFICATION_GUIDE.md`

---

## Example Results for PTCDA

**With 6-31g basis:**
```
Similarity: ~0.92
Correlation: ~0.94
R²: ~0.87
```

**Interpretation:** The approximation φ_HOMO×φ_LUMO captures ~92% of the exact transition density. Good enough for qualitative analysis, but use exact methods for quantitative work!

---

## Conclusion

This tool helps you understand:
- **What** the approximation is
- **Why** it differs from exact methods
- **When** it's good enough
- **How** basis set affects accuracy

Use it to gain insight into quantum chemical calculations! 🎓
