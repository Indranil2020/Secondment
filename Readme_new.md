# Transition Density Matrix Calculation with PySCF

This repository contains scripts for calculating transition density matrices and transition dipole moments of PTCDA using PySCF.

## Scripts

### 1. `tdm_calc_accurate.py` ⭐ **RECOMMENDED**
The most complete and accurate implementation:
- ✅ Full TDDFT (not TDA) for accurate transition dipoles
- ✅ Transition dipole moment calculations
- ✅ Oscillator strength calculations
- ✅ Natural Transition Orbital (NTO) analysis
- ✅ All three density types (transition, excited state, difference)
- ✅ Proper gauge origin handling
- ✅ Works with both test molecules and PTCDA

**Use this for production calculations and publication-quality results.**

### 2. `tdm_calc_test.py`
Basic implementation for testing:
- Uses TDA (Tamm-Dancoff Approximation)
- Density matrix calculations only
- Good for learning and quick tests

### 3. `tdm_calc_enhanced.py`
Enhanced density calculations:
- Multiple density visualization options
- Uses TDA
- No transition dipole calculations

## Quick Start

```bash
# Run the accurate calculation
python tdm_calc_accurate.py
```

For detailed information, see `CALCULATION_GUIDE.md`.

## Output Files

- `transition_density_state*.cube` - Transition density for dipole moments
- `excited_state_density_state*.cube` - Total excited state density
- `density_difference_state*.cube` - Density redistribution (best for visualization)
- `nto_state_*.molden` - Natural transition orbitals

## Visualization

**VMD:**
```bash
vmd density_difference_state1.cube
```

**Jmol:**
```bash
jmol density_difference_state1.cube
```

See `CALCULATION_GUIDE.md` for detailed visualization instructions.

## References

- PySCF Documentation: https://pyscf.org/user/tddft.html
- Example: https://github.com/pyscf/pyscf/blob/master/examples/1-advanced/030-transition_dipole.py
- Example: https://github.com/pyscf/pyscf/blob/master/examples/tddft/01-nto_analysis.py
