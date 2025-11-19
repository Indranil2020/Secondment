#!/usr/bin/env python
'''
Transition density matrix calculation using official PySCF methods.
Based on: https://github.com/pyscf/pyscf/blob/master/examples/tddft/22-density.py
'''

from pyscf import gto, dft, tdscf
from pyscf.tools import cubegen
import numpy as np

# 1. Define molecule (H2O)
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
td = tdscf.TDA(mf)
td.nstates = 5  # Number of excited states
td.kernel()

# 4. Calculate transition density using official PySCF method
state_id = 0  # First excited state (0-indexed)

def tda_density_matrix(td, state_id):
    '''
    Calculate excited state density matrix in AO basis.
    Official method from PySCF examples/tddft/22-density.py
    '''
    cis_t1 = td.xy[state_id][0]
    dm_oo = -np.einsum('ia,ka->ik', cis_t1.conj(), cis_t1)
    dm_vv = np.einsum('ia,ic->ac', cis_t1, cis_t1.conj())
    
    mf = td._scf
    dm = np.diag(mf.mo_occ)
    
    nocc = cis_t1.shape[0]
    dm[:nocc, :nocc] += dm_oo * 2
    dm[nocc:, nocc:] += dm_vv * 2
    
    mo = mf.mo_coeff
    dm = np.einsum('pi,ij,qj->pq', mo, dm, mo.conj())
    return dm

# Calculate excited state density
dm_excited = tda_density_matrix(td, state_id)

# Calculate density difference (excited - ground)
dm_diff = dm_excited - mf.make_rdm1()

# 5. Generate cube files
print(f"\nExcitation energy: {td.e[state_id]*27.211:.3f} eV")

# Excited state density
cubegen.density(mol, 'excited_state_density.cube', dm_excited, nx=80, ny=80, nz=80)
print("Generated: excited_state_density.cube")

# Density difference (recommended for visualization)
cubegen.density(mol, 'density_difference.cube', dm_diff, nx=80, ny=80, nz=80)
print("Generated: density_difference.cube")

print("\nVisualization tip:")
print("Use density_difference.cube in VMD with isovalues ±0.002")
print("Red (+) = electron gain, Blue (-) = electron loss")