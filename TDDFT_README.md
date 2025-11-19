# TDDFT Transition Density Matrix Calculation

## Overview
This directory contains scripts for calculating and visualizing transition density matrices from TDDFT (Time-Dependent Density Functional Theory) calculations using PySCF.

## Files

### Main Scripts
- **`tdm_calc_test.py`** - Simple script for basic transition density calculation
- **`tdm_calc_enhanced.py`** - Enhanced version that generates multiple types of density cube files

### Output Files
Three types of cube files are generated for each excited state:

1. **`excited_state_density_*.cube`** - Total electron density in the excited state
2. **`density_difference_*.cube`** - Difference between excited and ground state densities (RECOMMENDED)
3. **`transition_density_*.cube`** - Transition density matrix for transition moments

## Method

The code follows the **official PySCF approach** from:
- [PySCF examples/tddft/22-density.py](https://github.com/pyscf/pyscf/blob/master/examples/tddft/22-density.py)
- [PySCF User Guide: TDDFT](https://pyscf.org/user/tddft.html)

### Key Implementation Details

1. **TDA (Tamm-Dancoff Approximation)** is used instead of full TDDFT:
   ```python
   td = tdscf.TDA(mf)  # More stable than TDDFT
   ```

2. **Excited State Density Matrix** calculation:
   ```python
   def tda_density_matrix(td, state_id):
       cis_t1 = td.xy[state_id][0]
       dm_oo = -np.einsum('ia,ka->ik', cis_t1.conj(), cis_t1)
       dm_vv = np.einsum('ia,ic->ac', cis_t1, cis_t1.conj())
       
       dm = np.diag(mf.mo_occ)
       dm[:nocc, :nocc] += dm_oo * 2
       dm[nocc:, nocc:] += dm_vv * 2
       
       dm = np.einsum('pi,ij,qj->pq', mo, dm, mo.conj())
       return dm
   ```

3. **Density Difference** (excited - ground):
   ```python
   dm_diff = dm_excited - mf.make_rdm1()
   ```

## Usage

### Basic Usage
```bash
python tdm_calc_test.py
```

### Enhanced Version (Multiple States)
```bash
python tdm_calc_enhanced.py
```

## Visualization

### VMD (Recommended)

1. Load the density difference file:
   ```bash
   vmd density_difference_1.cube
   ```

2. In VMD GUI:
   - **Graphics → Representations**
   - Change **Drawing Method** to **Isosurface**
   - Create two representations:
     - **Rep 1**: Isovalue = `+0.002`, Color = **Red** (electron gain)
     - **Rep 2**: Isovalue = `-0.002`, Color = **Blue** (electron loss)

3. Adjust isovalue as needed (typical range: ±0.001 to ±0.005)

### Jmol

```tcl
isosurface ID "surf1" cutoff  0.002 density_difference_1.cube
isosurface ID "surf2" cutoff -0.002 density_difference_1.cube
```

## Understanding the Output

### Excited State Energies
The script prints excitation energies in eV:
```
State 1: 7.843 eV
State 2: 9.965 eV
State 3: 9.999 eV
```

### Density Interpretation
- **Red regions** (positive isosurface): Areas where electron density increases upon excitation
- **Blue regions** (negative isosurface): Areas where electron density decreases upon excitation
- This shows the **character of the electronic transition**

## Technical Notes

### Why TDA instead of TDDFT?
- TDA neglects the B matrix in the TDDFT equations
- More stable for systems with degeneracies or triplet instabilities
- Computationally less expensive
- Generally provides good accuracy for singlet excited states

### Normalization
- In TDA/RKS calculations, the amplitudes (X,Y) are normalized to 1/2
- Factor of 2 in the code accounts for spin-down contribution (closed-shell systems)

### Basis Set and Functional
Current settings:
- **Basis**: 6-31g (minimal basis for testing)
- **Functional**: B3LYP (hybrid functional)
- For production calculations, consider larger basis sets (e.g., 6-311+G(d,p), cc-pVTZ)

## References

1. PySCF Documentation: https://pyscf.org/
2. PySCF GitHub Examples: https://github.com/pyscf/pyscf/tree/master/examples/tddft
3. Martin, R. L., "Natural transition orbitals", J. Chem. Phys. 118, 4775-4777 (2003)

## Troubleshooting

### Issue: Only seeing H2O molecule, not transition density
**Solution**: This was the original problem. Make sure you're using the updated code that properly transforms the density matrix from MO to AO basis.

### Issue: Isosurface looks wrong in VMD
**Solution**: Adjust the isovalue. Start with ±0.002 and increase/decrease by factors of 2-5.

### Issue: Cube file is empty or has no data
**Solution**: Check that the TDDFT calculation converged. Look for "converged SCF energy" in the output.

## Example Workflow

```python
# 1. Define molecule
mol = gto.M(atom='O 0 0 0.117; H 0 0.757 -0.469; H 0 -0.757 -0.469', basis='6-31g')

# 2. Ground state DFT
mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.kernel()

# 3. TDDFT with TDA
td = tdscf.TDA(mf)
td.nstates = 5
td.kernel()

# 4. Calculate density
dm_excited = tda_density_matrix(td, state_id=0)
dm_diff = dm_excited - mf.make_rdm1()

# 5. Generate cube file
cubegen.density(mol, 'density_diff.cube', dm_diff, nx=80, ny=80, nz=80)
```

## Contact
For PySCF issues, see: https://github.com/pyscf/pyscf/issues
