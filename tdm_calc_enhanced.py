#!/usr/bin/env python
'''
Transition density matrix calculation for TDDFT excited states.
Based on PySCF official examples:
- examples/tddft/22-density.py
- User guide: https://pyscf.org/user/tddft.html
'''

from pyscf import gto, dft, tdscf
from pyscf.tools import cubegen
import numpy as np

# 1. Define molecule
mol = gto.M(
    atom = '''
    O  0.0000  0.0000  0.1173
    H  0.0000  0.7572 -0.4692
    H  0.0000 -0.7572 -0.4692
    ''',
    basis = '6-31g',
    verbose = 4
)

# 2. Ground state DFT calculation
mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.kernel()

# 3. TDDFT calculation using TDA (Tamm-Dancoff Approximation)
# TDA is more stable and commonly used for excited states
td = tdscf.TDA(mf)
td.nstates = 5
td.kernel()
td.analyze()  # Print excitation analysis

# Print all excitation energies
print("\n" + "="*60)
print("EXCITED STATE ENERGIES")
print("="*60)
for i, energy in enumerate(td.e):
    print(f"State {i+1}: {energy*27.211:.3f} eV")
print("="*60 + "\n")

# 4. Function to calculate transition density matrix
# Based on official PySCF example: examples/tddft/22-density.py
def tda_density_matrix(td, state_id):
    """
    Calculate the excited state density matrix in AO basis.
    Taking the TDA amplitudes as CIS coefficients.
    
    Parameters:
    -----------
    td : TDA or TDDFT object
        The TD calculation object
    state_id : int
        Excited state index (0-indexed)
    
    Returns:
    --------
    dm : ndarray
        Excited state density matrix in AO basis
    
    Reference:
    ----------
    PySCF examples/tddft/22-density.py
    """
    # Get TDA/CIS amplitudes for the excited state
    cis_t1 = td.xy[state_id][0]
    
    # Calculate density matrix changes in MO basis
    # dm_oo: occupied-occupied block (electron depletion)
    # dm_vv: virtual-virtual block (electron accumulation)
    dm_oo = -np.einsum('ia,ka->ik', cis_t1.conj(), cis_t1)
    dm_vv = np.einsum('ia,ic->ac', cis_t1, cis_t1.conj())
    
    # Start with ground state density matrix in MO basis
    mf = td._scf
    dm = np.diag(mf.mo_occ)
    
    # Add CIS contribution to density matrix
    nocc = cis_t1.shape[0]
    # Factor of 2 accounts for spin-down contribution (closed-shell)
    dm[:nocc, :nocc] += dm_oo * 2
    dm[nocc:, nocc:] += dm_vv * 2
    
    # Transform density matrix from MO basis to AO basis
    mo = mf.mo_coeff
    dm = np.einsum('pi,ij,qj->pq', mo, dm, mo.conj())
    
    return dm

def transition_density_matrix(td, state_id):
    """
    Calculate transition density matrix between ground and excited state.
    This represents the electronic transition, not the excited state density.
    
    Parameters:
    -----------
    td : TDA or TDDFT object
        The TD calculation object
    state_id : int
        Excited state index (0-indexed)
    
    Returns:
    --------
    dm_trans : ndarray
        Transition density matrix in AO basis
    """
    # Get TDA/CIS amplitudes
    cis_t1 = td.xy[state_id][0]
    
    # Construct transition density in MO basis
    # T_ia represents transition from occupied i to virtual a
    mf = td._scf
    nocc = cis_t1.shape[0]
    nmo = mf.mo_coeff.shape[1]
    
    dm_trans_mo = np.zeros((nmo, nmo))
    dm_trans_mo[:nocc, nocc:] = cis_t1
    dm_trans_mo[nocc:, :nocc] = cis_t1.T
    
    # Transform to AO basis
    mo = mf.mo_coeff
    dm_trans = np.einsum('pi,ij,qj->pq', mo, dm_trans_mo, mo.conj())
    
    return dm_trans

# 5. Generate cube files for selected states
states_to_visualize = [0, 1, 2]  # First three excited states (0-indexed)

print("\n" + "="*60)
print("GENERATING CUBE FILES")
print("="*60)

for state_id in states_to_visualize:
    print(f"\nState {state_id+1}: {td.e[state_id]*27.211:.3f} eV")
    
    # Option 1: Excited state density (total electron density in excited state)
    dm_excited = tda_density_matrix(td, state_id)
    filename_excited = f'excited_state_density_{state_id+1}.cube'
    cubegen.density(mol, filename_excited, dm_excited, nx=80, ny=80, nz=80)
    print(f"  Excited state density: {filename_excited}")
    
    # Option 2: Density difference (excited - ground state)
    dm_diff = dm_excited - mf.make_rdm1()
    filename_diff = f'density_difference_{state_id+1}.cube'
    cubegen.density(mol, filename_diff, dm_diff, nx=80, ny=80, nz=80)
    print(f"  Density difference: {filename_diff}")
    
    # Option 3: Transition density matrix (for transition moments)
    dm_trans = transition_density_matrix(td, state_id)
    filename_trans = f'transition_density_{state_id+1}.cube'
    cubegen.density(mol, filename_trans, dm_trans, nx=80, ny=80, nz=80)
    print(f"  Transition density: {filename_trans}")

print("\n" + "="*60)
print("VISUALIZATION INSTRUCTIONS")
print("="*60)
print("\nThree types of cube files generated:")
print("1. excited_state_density_*.cube - Total electron density in excited state")
print("2. density_difference_*.cube - Change in density (excited - ground)")
print("3. transition_density_*.cube - Transition density matrix")
print("\nRecommended: Use density_difference_*.cube for visualization")
print("\nIn VMD:")
print("  vmd density_difference_1.cube")
print("  Graphics > Representations > Drawing Method: Isosurface")
print("  Create two representations:")
print("    - Rep 1: Isovalue = +0.002, Color = Red (electron gain)")
print("    - Rep 2: Isovalue = -0.002, Color = Blue (electron loss)")
print("\nIn Jmol:")
print('  isosurface ID "surf1" cutoff  0.002 density_difference_1.cube')
print('  isosurface ID "surf2" cutoff -0.002 density_difference_1.cube')
print("="*60)
