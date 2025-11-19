# Script Comparison: Transition Density Matrix Calculations

## Feature Comparison Table

| Feature | tdm_calc_test.py | tdm_calc_enhanced.py | tdm_calc_accurate.py ⭐ |
|---------|------------------|----------------------|------------------------|
| **Method** | TDA | TDA | Full TDDFT |
| **Transition Dipole Moments** | ❌ No | ❌ No | ✅ Yes |
| **Oscillator Strengths** | ❌ No | ❌ No | ✅ Yes |
| **Gauge Origin Setting** | ❌ No | ❌ No | ✅ Yes |
| **NTO Analysis** | ❌ No | ❌ No | ✅ Yes |
| **Transition Density Matrix** | ❌ No | ✅ Yes | ✅ Yes |
| **Excited State Density** | ✅ Yes | ✅ Yes | ✅ Yes |
| **Density Difference** | ✅ Yes | ✅ Yes | ✅ Yes |
| **Cube File Generation** | ✅ Yes | ✅ Yes | ✅ Yes |
| **Molden File Output** | ❌ No | ❌ No | ✅ Yes |
| **Multiple States** | 1 state | 3 states | Configurable |
| **XYZ File Support** | ❌ No | ❌ No | ✅ Yes |
| **Documentation** | Basic | Good | Comprehensive |
| **Publication Ready** | ❌ No | ❌ No | ✅ Yes |

## When to Use Each Script

### Use `tdm_calc_test.py` when:
- Learning PySCF basics
- Quick density visualization tests
- Don't need transition dipoles
- Working with simple molecules

### Use `tdm_calc_enhanced.py` when:
- Need multiple density types
- Want detailed density analysis
- Don't need transition dipoles or NTO
- Prefer TDA for stability

### Use `tdm_calc_accurate.py` when: ⭐
- Need accurate transition dipole moments
- Calculating oscillator strengths
- Performing NTO analysis
- Publishing results
- Working with PTCDA or large molecules
- Need complete excited state characterization

## Key Technical Differences

### 1. TDDFT Method

**TDA (Tamm-Dancoff Approximation):**
```python
td = tdscf.TDA(mf)  # Used in test and enhanced scripts
```
- Only X amplitudes (Y = 0)
- Faster computation
- Less accurate for transition dipoles
- Good for density visualization

**Full TDDFT:**
```python
td = tddft.TDDFT(mf)  # Used in accurate script
```
- Both X and Y amplitudes
- More accurate transition properties
- Required for proper oscillator strengths
- Recommended for spectroscopy

### 2. Transition Density Matrix Construction

**TDA version (test/enhanced):**
```python
# Only uses X amplitudes
cis_t1 = td.xy[state_id][0]
dm_trans_mo[:nocc, nocc:] = cis_t1
```

**Full TDDFT version (accurate):**
```python
# Uses both X and Y amplitudes
X, Y = td.xy[state_id]
dm_trans_mo[:nocc, nocc:] = (X + Y).reshape(nocc, nvir)
```

### 3. Transition Dipole Calculation

**Not implemented in test/enhanced scripts.**

**Accurate script:**
```python
def calculate_transition_dipole(td, state_id):
    X, Y = td.xy[state_id]
    # Construct transition density matrix
    t_dm1_ao = construct_transition_density(X, Y)
    # Calculate dipole: μ = Tr(μ_op * T)
    tdm = np.einsum('xij,ji->x', dip_ints, t_dm1_ao)
    return tdm
```

### 4. Gauge Origin

**Test/Enhanced scripts:** No gauge origin setting (origin at [0,0,0])

**Accurate script:** Sets gauge origin to nuclear charge center
```python
nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
mol.set_common_orig_(nuc_charge_center)
```

This is crucial for meaningful transition dipole moments.

### 5. NTO Analysis

**Not in test/enhanced scripts.**

**Accurate script:**
```python
weights, nto_coeff = td.get_nto(state=i+1, verbose=4)
molden.from_mo(mol, f'nto_state_{i+1}.molden', nto_coeff)
```

NTOs provide intuitive visualization of hole and particle orbitals.

## Output Comparison

### tdm_calc_test.py
```
Output:
- excited_state_density.cube (1 state)
- density_difference.cube (1 state)
```

### tdm_calc_enhanced.py
```
Output:
- excited_state_density_1.cube (3 states)
- density_difference_1.cube (3 states)
- transition_density_1.cube (3 states)
```

### tdm_calc_accurate.py ⭐
```
Output:
- transition_density_state1.cube (configurable states)
- excited_state_density_state1.cube
- density_difference_state1.cube
- nto_state_1.molden
- Transition dipole moments (printed)
- Oscillator strengths (printed)
- NTO weights (printed)
```

## Performance Comparison

For PTCDA (38 atoms, 6-31g basis):

| Script | Ground State | Excited States | Total Time |
|--------|--------------|----------------|------------|
| test (TDA) | ~2 min | ~5 min | ~7 min |
| enhanced (TDA) | ~2 min | ~5 min | ~7 min |
| accurate (TDDFT) | ~2 min | ~15 min | ~17 min |

**Note:** Full TDDFT is ~3x slower than TDA but provides accurate transition properties.

## Accuracy Comparison

### Transition Dipole Moments

| Method | Accuracy | Use Case |
|--------|----------|----------|
| TDA (test/enhanced) | Not calculated | N/A |
| Full TDDFT (accurate) | High accuracy | Spectroscopy, comparison with experiment |

### Oscillator Strengths

| Method | Accuracy | Use Case |
|--------|----------|----------|
| TDA | Approximate | Qualitative trends |
| Full TDDFT | Accurate | Quantitative comparison with experiment |

### Excitation Energies

| Method | Typical Error | Use Case |
|--------|---------------|----------|
| TDA | 0.1-0.3 eV | Good for most cases |
| Full TDDFT | 0.1-0.2 eV | Slightly better |

## Migration Guide

### From test.py to accurate.py

1. Change import:
```python
# Old
from pyscf import tdscf
td = tdscf.TDA(mf)

# New
from pyscf import tddft
td = tddft.TDDFT(mf)
```

2. Update density matrix function to use X and Y:
```python
# Old
cis_t1 = td.xy[state_id][0]

# New
X, Y = td.xy[state_id]
```

3. Add transition dipole calculation (see accurate script)

### From enhanced.py to accurate.py

Same as above, plus:
- Add NTO analysis
- Add oscillator strength calculation
- Set gauge origin

## Recommendations

### For Learning
Start with `tdm_calc_test.py` to understand basics.

### For Research
Use `tdm_calc_accurate.py` for all production calculations.

### For Quick Checks
Use `tdm_calc_enhanced.py` if you only need density visualization and TDA is sufficient.

### For Publications
**Always use `tdm_calc_accurate.py`** to ensure:
- Accurate transition dipole moments
- Proper oscillator strengths
- Complete excited state characterization
- Reproducible results

## Validation Checklist

When using any script, verify:

- [ ] Excitation energies are reasonable (compare with literature)
- [ ] Oscillator strengths match experimental absorption
- [ ] NTO weights indicate single-excitation character (>0.5)
- [ ] Density difference integrates to ~0
- [ ] Transition dipoles are gauge-origin independent (for accurate script)

## References

1. **PySCF Examples:**
   - `examples/tddft/22-density.py` - Density matrices
   - `examples/1-advanced/030-transition_dipole.py` - Transition dipoles
   - `examples/tddft/01-nto_analysis.py` - NTO analysis

2. **Theory:**
   - Casida, M. E. "Time-dependent density functional response theory"
   - Hirata & Head-Gordon, "TDDFT within the Tamm–Dancoff approximation"

3. **PySCF Documentation:**
   - https://pyscf.org/user/tddft.html
   - https://pyscf.org/pyscf_api_docs/pyscf.tddft.html
